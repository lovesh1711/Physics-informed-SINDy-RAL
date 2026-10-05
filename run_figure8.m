%% ========================================================================
%  run_figure8.m
%  Closed-loop figure-8 tracking with THREE controllers, to demonstrate that
%  the beta_z diagnostic predicts the Koopman failure MODE:
%    K-narrow : Koopman trained on hover data  -> beta_z < 0 (correct sign,
%               frozen-upright direction)
%    K-broad  : Koopman trained on broad data  -> beta_z > 0 (inverted sign)
%    SINDy    : physics-informed state-space model (the cure)
%  Reference: differentially-flat figure-8 (gen_figure8_ref.m).
% ========================================================================
clc; clear; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'utils'),   '-begin');
addpath(fullfile(thisDir,'koopman'), '-begin');
addpath(fullfile(thisDir,'sindy'),   '-begin');

params.m = 4.34; params.J = diag([0.0820;0.0845;0.1377]); params.g = 9.81;
u_min = [10;-30;-30;-30]; u_max = [300;30;30;30];
ts = 0.01; dt = 0.001; steps = round(ts/dt);
n = 13; nu = 4; mg = params.m*params.g; Nh = 7;
p_lift = 3; Nz = n + 9*p_lift; C_koop = [eye(n), zeros(n,9*p_lift)];

%% ---------- figure-8 reference ----------
T_total = 8; A = 2.0; w = 1.6; z0 = 0;
[Xref, Uref, tvec] = gen_figure8_ref(params, ts, T_total, A, w, z0);
Nsim = size(Xref,2) - Nh - 1;                 % leave room for the horizon
fprintf('Figure-8: %d control steps, thrust %.1f-%.1f N, max bank %.0f deg\n', ...
        Nsim, min(Uref(1,:)), max(Uref(1,:)), fig8_maxbank(Xref));

%% ---------- train the two Koopman models ----------
[Xn,Xpn,Un,Xdn] = gen_train(params,ts,dt,steps,mg,u_min,u_max,'hover',1);
[A_kn,B_kn] = train_edmd(Xn,Xpn,Un,p_lift,n);
bz_n = B_kn(6,1);
[Xb,Xpb,Ub] = gen_train(params,ts,dt,steps,mg,u_min,u_max,'broad',11);
[A_kb,B_kb] = train_edmd(Xb,Xpb,Ub,p_lift,n);
bz_b = B_kb(6,1);

%% ---------- train SINDy (hover data, finite-difference derivatives) ----------
[Theta,~] = sindy_library(Xn,Un);
Xi = stlsq(Theta, Xdn, 0.10, 10);           % same settings as run_comparison.m
fprintf('beta_z : K-narrow = %+.3e (%.2f e-3),  K-broad = %+.3e (%.2f e-3)\n', ...
        bz_n,1e3*bz_n, bz_b,1e3*bz_b);
fprintf('SINDy: nnz(Xi)=%d\n', nnz(Xi));

%% ---------- closed-loop runs ----------
Xk_n = run_koopman(A_kn,B_kn,C_koop,Nz,p_lift,Xref,Uref,params,ts,dt,steps,Nh,n,nu,Nsim);
Xk_b = run_koopman(A_kb,B_kb,C_koop,Nz,p_lift,Xref,Uref,params,ts,dt,steps,Nh,n,nu,Nsim);
Xs   = run_sindy(Xi,Xref,Uref,params,ts,dt,steps,Nh,n,nu,Nsim,u_min,u_max);

%% ---------- metrics ----------
% logged state k is at t = k*ts, i.e. reference sample k+1
rmse = @(X) sqrt(mean(sum((X(1:3,1:Nsim)-Xref(1:3,2:Nsim+1)).^2,1)));
fprintf('\n--- position tracking RMSE (m) ---\n');
fprintf('  K-narrow (beta_z<0): %.3f\n', rmse(Xk_n));
fprintf('  K-broad  (beta_z>0): %.3f\n', rmse(Xk_b));
fprintf('  SINDy             : %.3f\n', rmse(Xs));

%% ---------- save + quick plot ----------
resultsDir = fullfile(thisDir,'results');
save(fullfile(resultsDir,'figure8_results.mat'), ...
     'Xref','Uref','tvec','Xk_n','Xk_b','Xs','bz_n','bz_b','Nsim','-v7.3');

figure('Color','w','Position',[80 80 560 520]);
plot(Xref(2,1:Nsim),Xref(1,1:Nsim),'k-','LineWidth',2); hold on;
plot(Xk_n(2,1:Nsim),Xk_n(1,1:Nsim),'-','LineWidth',1.3);
plot(Xk_b(2,1:Nsim),Xk_b(1,1:Nsim),'-','LineWidth',1.3);
plot(Xs(2,1:Nsim), Xs(1,1:Nsim), '-','LineWidth',1.3);
axis equal; grid on; xlabel('East y [m]'); ylabel('North x [m]');
legend('reference','K-narrow (\beta_z<0)','K-broad (\beta_z>0)','SINDy','Location','best');
title('Figure-8 tracking');
exportgraphics(gcf, fullfile(resultsDir,'figure8_quicklook.png'),'Resolution',200);
fprintf('saved results/figure8_quicklook.png\nDONE\n');

%% ===================== local functions ==============================
function b = fig8_maxbank(Xref)
    e3=[0;0;1]; b=0;
    for k=1:size(Xref,2)
        R=QuatToRot(Xref(7:10,k)); v=R'*e3; b=max(b,acosd(min(max(v(3),-1),1)));
    end
end

function [X,Xp,U,Xdot] = gen_train(params,ts,dt,steps,mg,u_min,u_max,mode,seed)
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));
    n=13; nu=4;
    if strcmp(mode,'hover')
        Ntraj=100; T=0.5; Sig=diag([10;10;10;10]); hold_len=5;
    else
        Ntraj=600; T=1.0; Sig=diag([15;15;15;15]); hold_len=-1;
    end
    Ns=round(T/ts); X=zeros(n,Ntraj*Ns); Xp=X; Xdot=X; U=zeros(nu,Ntraj*Ns); kk=0;
    for tr=1:Ntraj
        if strcmp(mode,'hover')
            x=[zeros(6,1);1;0;0;0;zeros(3,1)];
        else
            p0=1.0*(2*rand(3,1)-1); v0=2.0*(2*rand(3,1)-1);
            ax=randn(3,1); ax=ax/max(norm(ax),1e-9); ang=(pi/4)*rand;
            q0=normalize_quat([cos(ang/2);sin(ang/2)*ax]); w0=3.0*(2*rand(3,1)-1);
            x=[p0;v0;q0;w0];
        end
        x(7:10)=normalize_quat(x(7:10));
        if hold_len==-1, if mod(tr,2)==0, hl=3; else, hl=7; end; else, hl=hold_len; end
        u=[mg;0;0;0];
        for k=1:Ns
            if (k==1)||(mod(k-1,hl)==0)
                u=mvnrnd(zeros(4,1),Sig).'; u(1)=u(1)+mg; u=min(max(u,u_min),u_max);
            end
            xk=x; xs=zeros(n,2);
            for j=1:steps
                x=rk4_step(x,u,params,dt); x(7:10)=normalize_quat(x(7:10));
                if j<=2, xs(:,j)=x; end
            end
            kk=kk+1; X(:,kk)=xk; Xp(:,kk)=x; U(:,kk)=u;
            for j=1:2, if dot(xs(7:10,j),xk(7:10))<0, xs(7:10,j)=-xs(7:10,j); end, end  % align q sign
            Xdot(:,kk)=(-3*xk+4*xs(:,1)-xs(:,2))/(2*dt);   % 2nd-order forward FD
        end
    end
end

function [A_koop,B_koop] = train_edmd(X,Xp,U,p_lift,n)
    K=size(X,2); Nz=n+9*p_lift; Z=zeros(Nz,K); Zp=zeros(Nz,K);
    for k=1:K
        xk=X(:,k); xk(7:10)=normalize_quat(xk(7:10));
        xp=Xp(:,k); xp(7:10)=normalize_quat(xp(7:10));
        Z(:,k)=observable_phi(xk,p_lift); Zp(:,k)=observable_phi(xp,p_lift);
    end
    ZU=[Z;U]; G=Zp*ZU'/(ZU*ZU'+1e-6*eye(size(ZU,1)));
    A_koop=G(:,1:Nz); B_koop=G(:,Nz+1:end);
end

function Xcl = run_koopman(A_koop,B_koop,C_koop,Nz,p_lift,Xref,Uref,params,ts,dt,steps,Nh,n,nu,Nsim)
    % Koopman TV-MPC with integral action (mirrors run_comparison.m)
    Q_pos=diag([500 500 500]); Q_vel=diag([50 50 50]);
    Q_q=diag([10 10 10 10]);   Q_w=diag([100 100 100]);
    S_int=1*eye(3); eta_max=0.5; Nz_aug=Nz+3+1; C_p=[eye(3),zeros(3,n-3)]*C_koop;
    Qtk=zeros(Nz,Nz); Qtk(1:13,1:13)=blkdiag(Q_pos,Q_vel,Q_q,Q_w);
    Qta=blkdiag(Qtk,S_int,0); Qbar=kron(eye(Nh),Qta);
    Ru=diag([0.01,5,5,5]); Rbar=kron(eye(Nh),Ru);
    umn=[10;-30;-30;-30]; umx=[300;30;30;30];   % same bounds as the SINDy MPC
    Ulb=repmat(umn,Nh,1); Uub=repmat(umx,Nh,1);
    opts=optimoptions('quadprog','Display','none','Algorithm','interior-point-convex');
    Nref=size(Xref,2);
    zref=zeros(Nz,Nref); for k=1:Nref, zref(:,k)=observable_phi([Xref(1:6,k);normalize_quat(Xref(7:10,k));Xref(11:13,k)],p_lift); end
    pref=Xref(1:3,:);
    x=[Xref(1:6,1);normalize_quat(Xref(7:10,1));Xref(11:13,1)];
    z=observable_phi(x,p_lift); eta=zeros(3,1); zaug=[z;eta;1];
    Xcl=zeros(n,Nsim);
    for k=1:Nsim
        rk=min(k,Nref);
        Aseq=cell(Nh,1); Bseq=cell(Nh,1);
        for i=1:Nh
            idx=min(rk+i-1,Nref);
            [Ai,Bi]=local_koopman_linearization(zref(:,idx),Uref(:,idx),A_koop,B_koop);
            prf=pref(:,idx);
            Aseq{i}=[Ai, zeros(Nz,3), zeros(Nz,1); ts*C_p, eye(3), -ts*prf; zeros(1,Nz), zeros(1,3), 1];
            Bseq{i}=[Bi; zeros(3,nu); zeros(1,nu)];
        end
        Aqp=zeros(Nz_aug*Nh,Nz_aug); Bqp=zeros(Nz_aug*Nh,nu*Nh); Ap=eye(Nz_aug);
        for i=1:Nh
            Ap=Aseq{i}*Ap; Aqp((i-1)*Nz_aug+1:i*Nz_aug,:)=Ap;
            for j=1:i
                Apr=eye(Nz_aug); for s=j+1:i, Apr=Aseq{s}*Apr; end
                Bqp((i-1)*Nz_aug+1:i*Nz_aug,(j-1)*nu+1:j*nu)=Apr*Bseq{j};
            end
        end
        Y=zeros(Nz_aug*Nh,1);
        for i=1:Nh
            idx=min(rk+i-1,Nref); b=(i-1)*Nz_aug;
            Y(b+1:b+Nz)=zref(:,idx); Y(b+Nz+4)=1;
        end
        H=2*(Bqp'*Qbar*Bqp+Rbar); H=0.5*(H+H')+1e-8*eye(size(H));
        f=2*(Bqp'*Qbar*(Aqp*zaug-Y));
        [U0,~,ef]=quadprog(H,f,[],[],[],[],Ulb,Uub,[],opts);
        if ef<=0||any(isnan(U0)), U0=min(max(-(H\f),Ulb),Uub); end
        u0=min(max(U0(1:nu),umn),umx);
        for j=1:steps, x=rk4_step(x,u0,params,dt); x(7:10)=normalize_quat(x(7:10)); end
        z=observable_phi(x,p_lift);
        eta=max(min(eta+(x(1:3)-pref(:,min(rk+1,Nref)))*ts,eta_max),-eta_max);   % same time t_{k+1}
        zaug=[z;eta;1]; Xcl(:,k)=x;
    end
end

function Xcl = run_sindy(Xi,Xref,Uref,params,ts,dt,steps,Nh,n,nu,Nsim,u_min,u_max)
    Qx=blkdiag(diag([300 200 1]),diag([100 100 100]),eye(4),diag([2500 2500 2500]));
    Qbar=kron(eye(Nh),Qx); Rbar=kron(eye(Nh),diag([0.2,4,4,4]));
    opts=optimoptions('quadprog','Display','none','Algorithm','interior-point-convex');
    Nref=size(Xref,2);
    x=[Xref(1:6,1);normalize_quat(Xref(7:10,1));Xref(11:13,1)];
    Xcl=zeros(n,Nsim);
    for k=1:Nsim
        rk=min(k,Nref); dx0=x-Xref(:,rk);
        Aseq=cell(Nh,1); Bseq=cell(Nh,1); DUlb=zeros(nu*Nh,1); DUub=zeros(nu*Nh,1);
        for i=1:Nh
            idx=min(rk+i-1,Nref);
            [Ac,Bc]=sindy_jacobian_fd(Xref(:,idx),Uref(:,idx),Xi);
            Aseq{i}=eye(n)+ts*Ac; Bseq{i}=ts*Bc;
            b=(i-1)*nu; DUlb(b+1:b+nu)=u_min-Uref(:,idx); DUub(b+1:b+nu)=u_max-Uref(:,idx);
        end
        Aqp=zeros(n*Nh,n); Bqp=zeros(n*Nh,nu*Nh); Ap=eye(n);
        for i=1:Nh
            Ap=Aseq{i}*Ap; Aqp((i-1)*n+1:i*n,:)=Ap;
            for j=1:i
                Apr=eye(n); for s=j+1:i, Apr=Aseq{s}*Apr; end
                Bqp((i-1)*n+1:i*n,(j-1)*nu+1:j*nu)=Apr*Bseq{j};
            end
        end
        H=2*(Bqp'*Qbar*Bqp+Rbar); H=0.5*(H+H')+1e-9*eye(size(H));
        f=2*(Bqp'*Qbar*(Aqp*dx0));
        [DU,~,ef]=quadprog(H,f,[],[],[],[],DUlb,DUub,[],opts);
        if ef<=0||any(isnan(DU)), DU=min(max(-(H\f),DUlb),DUub); end
        u0=min(max(Uref(:,rk)+DU(1:nu),u_min),u_max);
        for j=1:steps, x=rk4_step(x,u0,params,dt); x(7:10)=normalize_quat(x(7:10)); end
        Xcl(:,k)=x;
    end
end
