function [sim,mpc_data] = mpc_simulation(G,DG,xDG0,xe0,dsim,param,opts)
tic
zmean = param.zmean;
mu = param.mu;
nc = length(mu);
muo = param.muo;
rho = param.rho;
gam = param.gam;
gamo = param.gamo;
dmax = param.dmax;
cp0 = param.cp0;
cf0 = param.cf0;
Q = param.Q;
N = param.N;
K = param.K;
L = param.L;
alpha = rho.^2./(1-rho.^2);
beta = alpha.*(gam-mu);
betao = alpha.*(gamo-muo);
pG = G.InputGroup.p;
dG = G.InputGroup.d;
uG = G.InputGroup.u;
qG = G.OutputGroup.q;
zG = G.OutputGroup.z;
yG = G.OutputGroup.y;
yyG = G.OutputGroup.y;
nx = length(G.A);
GK = lft(G([qG zG yG yyG],[pG dG uG uG]),K);%(:,yK));
nGK = length(GK.A);
nu = length(uG);
nq = length(qG);
nz = length(zG);
%ny = length(yG);
ny = length(yyG);
uGK = GK.InputGroup.u;
Dscript = orth(GK.D(qG,uGK));
Cscript = (eye(nq) - Dscript*Dscript')*GK.C(qG,:); % Note that orth yields 
% an ortonormal basis, i.e., Dscript'*Dscript = eye and hence not needed.

% Get terminal ingredients:
[S,T,cff,GKxnff] = terminal_conditions(Q,GK.A,GK.C(qG,:),GK.C(zG,:),param);

%% Definition of MPC OCP
opti = casadi.Opti();
if param.only_K
    init_opt = 1;
else
    init_opt = opti.variable(1);  % parameter that optimizes the initial condition
    opti.subject_to(init_opt>=0);
    opti.subject_to(init_opt<=1);
end
unf = opti.variable(nu,N);       % nominal input trajectory
xGKnf_pred = opti.parameter(nGK,1);  % predicted state
xGKnf_obs = opti.parameter(nGK,1);   % observed state
cf_pred = opti.parameter(nc,1);  % predicted tube size
cf_obs = opti.parameter(nc,1);   % observed tube size

xGKnf = casadi.MX(zeros(nGK,N+1));
qnf = casadi.MX(zeros(nq,N));
znf = casadi.MX(zeros(nz,N));
cf = opti.variable(nc,N+1);
xGKnf(:,1)=xGKnf_pred*init_opt+(1-init_opt)*xGKnf_obs;
opti.subject_to(cf(:,1)==cf_pred*init_opt+(1-init_opt)*cf_obs);
ell = casadi.MX(0);
for k=1:N
    xGKnf(:,k+1) = GK.A*xGKnf(:,k)+GK.B(:,uGK)*unf(:,k);
    qnf(:,k) = GK.C(qG,:)*xGKnf(:,k)+GK.D(qG,uGK)*unf(:,k);
    znf(:,k) = GK.C(zG,:)*xGKnf(:,k)+GK.D(zG,uGK)*unf(:,k);
    qdmax2 = (qnf(:,k)'*qnf(:,k)+dmax^2);
    opti.subject_to(cf(:,k+1) == rho.^2.*cf(:,k) + mu.*rho.^2*qdmax2);
    const_tight = gam./alpha.*(cf(:,k) + beta*qdmax2);
    opti.subject_to(znf(:,k)-zmean<=1);
    opti.subject_to(znf(:,k)-zmean>=-1);
    opti.subject_to(const_tight<=(1+zmean-znf(:,k)).^2);
    opti.subject_to(const_tight<=(znf(:,k)-zmean+1).^2);
    ell = ell + [xGKnf(:,k)' unf(:,k)']*Q*[xGKnf(:,k); unf(:,k)];
end
ell = ell + xGKnf(:,k)'*S*xGKnf(:,k);
opti.minimize(ell);
% terminal conditions:
TT = T'*T;
figure(5)
E = ellipsoid([0;0], inv(TT(1:2,1:2))*GKxnff);
plot(E,'b');
opti.subject_to(xGKnf(:,end)'*TT*xGKnf(:,end) <= GKxnff )
opti.subject_to(cf(:,end) <= cff )
p_opts = struct('expand', true);
s_opts = struct('max_iter',1e4, 'print_level', 0);
opti.solver('ipopt', p_opts, s_opts);
time_setup = toc;
%% Simulation
tic
Tsim = opts.Tsim;
%DG = lft(Delta,G); % real sys
nDG = length(DG.A);
dDG = DG.InputGroup.d;
uDG = DG.InputGroup.u;
xDGsim = zeros(nDG,Tsim+1);
nK = length(K.A);
xKsim = zeros(nK,Tsim+1);
nL = length(L.A);
xLsim = zeros(nL,Tsim+1);
xGKesim = zeros(nx+nK,Tsim);
xGKnpsim = zeros(nx+nK,Tsim+1);
xGKnpsim(1:nx,1) = xe0;
mpc_data = struct();
usim = zeros(nu,Tsim);
eunfsim = zeros(nu,Tsim);
zsim = zeros(nz,Tsim);
qnpsim = zeros(nq,Tsim);
ynpsim = zeros(ny,Tsim);
cpsim = zeros(nc,Tsim);
cpsim(:,1) = cp0;
zDG = DG.OutputGroup.z;
yyDG = DG.OutputGroup.y;
ysim = zeros(ny,Tsim);
xDGsim(:,1) = xDG0;
xGKe0 = [xe0; zeros(nK,1)];
opti.set_value(xGKnf_pred, xGKe0)
opti.set_value(xGKnf_obs, xGKe0)
for t=1:Tsim
    % measure
    ysim(:,t) = DG.C(yyDG,:)*xDGsim(:,t) + DG.D(yyDG,dDG)*dsim(:,t);
    % observe
    ynpsim(:,t) = GK.C(yyG,:)*xGKnpsim(:,t);
    eynp = ysim(:,t)-ynpsim(:,t); 
    xLsim(:,t+1) = L.A*xLsim(:,t) + L.B*eynp;
    exGKnpe = L.C*xLsim(:,t) + L.D*eynp;
    xGKesim(:,t) = xGKnpsim(:,t) + exGKnpe;
    % prepare OCP
    opti.set_value( xGKnf_obs, xGKesim(:,t) )
    Cth = Cscript*xGKnpsim(:,t);
    opti.set_value( cf_obs, cpsim(:,t)+betao*(Cth'*Cth+dmax^2) )
    if t==1
        opti.set_value( xGKnf_pred, xGKe0 )
        opti.set_value( cf_pred, cf0 )
    else
        opti.set_value( xGKnf_pred, mpc_data(t-1).xGKnf(:,2) )
        opti.set_value( cf_pred, mpc_data(t-1).cf(:,2) )
        opti.set_initial( init_opt, 1 )
        opti.set_initial( unf, [mpc_data(t-1).unf(:,2:end) zeros(nu,1)])
    end
    try
        % solve OCP
        sol=opti.solve;
        % get solution of OCP
        mpc_data(t).unf = sol.value(unf);
        mpc_data(t).init_opt = sol.value(init_opt);
        mpc_data(t).xGKnf = sol.value(xGKnf);
        mpc_data(t).znf = sol.value(znf);
        mpc_data(t).qnf = sol.value(qnf);
        mpc_data(t).cf = sol.value(cf);
    catch
        if t==1
            error("Problem seems to be infeasible at t=0")
        end
        disp('OCP not solved, continue with warm-sart')
        mpc_data(t).unf = [mpc_data(t-1).unf(:,2:end) zeros(nu,1)];
        mpc_data(t).init_opt = 1;
        mpc_data(t).xGKnf = [mpc_data(t-1).xGKnf(:,2:end) GK.A*mpc_data(t-1).xGKnf(:,end)];
        mpc_data(t).znf = [mpc_data(t-1).znf(:,2:end) GK.C(zG,:)*mpc_data(t-1).xGKnf(:,end)];
        mpc_data(t).qnf = [mpc_data(t-1).qnf(:,2:end) GK.C(qG,:)*mpc_data(t-1).xGKnf(:,end)];
        qdmax2 = mpc_data(t).qnf(:,end)'*mpc_data(t).qnf(:,end)+dmax^2;
        mpc_data(t).cf = [mpc_data(t-1).cf(:,2:end) rho.^2.*mpc_data(t-1).cf(:,2:end)+mu.*rho.^2*qdmax2];
    end
    % compute actual input
    xKsim(:,t+1) = K.A*xKsim(:,t) + K.B*ysim(:,t);
    eunfsim(:,t) = K.C*xKsim(:,t) + K.D*ysim(:,t);
    usim(:,t) = mpc_data(t).unf(:,1) + eunfsim(:,t);
    % apply input to system
    xDGsim(:,t+1) = DG.A*xDGsim(:,t) + DG.B(:,dDG)*dsim(:,t) + DG.B(:,uDG)*usim(:,t);
    zsim(:,t) = DG.C(zDG,:)*xDGsim(:,t) + DG.D(zDG,[dDG uDG])*[dsim(:,t); usim(:,t)];
    % apply input to xGKnp
    xGKnpsim(:,t+1) = GK.A*xGKnpsim(:,t) + GK.B(:,uGK)*mpc_data(t).unf(:,1);
    qnpsim(:,t) = GK.C(qG,:)*xGKnpsim(:,t) + GK.D(qG,uGK)*mpc_data(t).unf(:,1);
    cpsim(:,t+1) = rho.^2.*cpsim(:,t) + muo.*rho.^2*(qnpsim(:,t)'*qnpsim(:,t)+dmax^2);
    xLsim(:,t+1) = L.A*xLsim(:,t) + L.B*eynp;
    exGKnpe = L.C*xLsim(:,t) + L.D*eynp;
    xGKesim(:,t) = xGKnpsim(:,t) + exGKnpe;
end
time_sim = toc;
sim = struct();
sim.xDG = xDGsim;
sim.xGKnp = xGKnpsim;
sim.qnp = qnpsim;
sim.xGKe = xGKesim;
sim.cp = cpsim;
sim.xL = xLsim;
sim.xK = xKsim;
sim.u = usim;
sim.z = zsim;
sim.eunf = eunfsim;
sim.xGKnf = zeros(nGK,Tsim);
sim.unf = zeros(nu,Tsim);
sim.znf = zeros(nz,Tsim);
sim.qnf = zeros(nq,Tsim);
sim.cf = zeros(nc,Tsim);
sim.init_opt = zeros(1,Tsim);
for t = 1:Tsim
    sim.xGKnf(:,t) = mpc_data(t).xGKnf(:,1);
    sim.unf(:,t) = mpc_data(t).unf(:,1);
    sim.init_opt(t) = mpc_data(t).init_opt(1);
    sim.znf(:,t) = mpc_data(t).znf(:,1);
    sim.qnf(:,t) = mpc_data(t).qnf(:,1);
    sim.cf(:,t) = mpc_data(t).cf(:,1);
end
sim.time_sim = time_sim;
sim.time_setup = time_setup;

