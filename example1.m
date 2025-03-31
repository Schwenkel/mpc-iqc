clear all
addpath("multi-objective-iqc-synthesis")
%% System
tic
% define model
m = 1;  % mass
c = 1;  % spring constant
d = 1;  % damper constant
Ts=0.1; % sampling time
Ac = [0 1; -c/m -d/m];   % continuous-time dynamics
nx=length(Ac);
Bc = [0 0 0; 1 1/m 1/m]; % continuous-time dynamics
AB = expm([Ac Bc; zeros(3,5)]*Ts);
Ad = AB(1:2,1:2);   % discrete-time dynamics
Bd = AB(1:2,3:5);   % discrete-time dynamics
C = [ c/m d/m;  % q
         1 0;   % z1 = x1
         0 1;   % z2 = x2
         0 0;   % z3 = u
         1 0;   % y1 = x1
         0 1];  % y2 = x2
p = 1;   np=length(p);  % input p
d = 2;   nd=length(d);  % input d
u = 3;   nu=length(u);  % input u
q = 1;   nq=length(q);  % output q
z = 2:4; nz=length(z);  % output z
y = 5:6; ny=length(y);  % output y
D = zeros(y(end),u(end));
D(q,u) = -1/m;  % q = -1/m*u+...
D(z(3),u)=1;    % z3=u;
% constraints
zmin = [-0.1; -0.25;-1];
zmax = [1; 0.05;1];
zscaling = (zmax-zmin)/2;    
zmean = (zmax+zmin)/2./zscaling; 
% then z./zscaling-zmean must be in [-1,1]
% Rescale:
dscaling = 0.3;  % a means to trade off the influence of d and q on the 
                 % peak-to-peak gain from w=[d; q] to z.
Bd(:,d) = Bd(:,d)*dscaling;
D(:,d) = D(:,d)*dscaling;
C(z,:) = diag(1./zscaling)*C(z,:);
D(z,:) = diag(1./zscaling)*D(z,:);
% use w = [ d; q ] as input
Bd_w= [ Bd(:,p) zeros(nx,nq) Bd(:,d) Bd(:,u) ];
D_w = [ D(q,p)        eye(nq)       D(q,[d u]);
        D(z,p)        zeros(nz,nq)  D(z,[d u]);
        D(y,p)        zeros(ny,nq)  D(y,[d u])];
w = d(1):(d(end)+nq);
G_w = ss(Ad,Bd_w,C,D_w,Ts);
%G_K = ss(Ad,Bd,C,D,Ts);
%mdm = m_delta_max;
%Wq = -mdm*4*(zs-1)^2/(mdm*4*(zs-1)^2+2*Ts*d2*(zs^2-1)+Ts^2*c2*(zs+1)^2);
%G_K = blkdiag(Wq,eye(nz+ny+nq))*G_K; 
G_w.InputGroup.p = p;
G_w.InputGroup.w = w;
G_w.InputGroup.d = d+nq;
G_w.InputGroup.u = u+nq;
G_w.OutputGroup.q = q;
G_w.OutputGroup.z = z;
G_w.OutputGroup.y = y;
%G = balreal(G);
%% Uncertainty and IQC:
mdm = 0.1;         % True value m_delta in [-mdm, mdm]
c2 = 1; d2 = 0.5;  % spring and damper constant of uncertain block
zs=tf('z',Ts);
% worst case Delta discretized using Tustin method:
Delta_max = (c2*Ts^2*(zs+1)^2+d2*2*(zs-1)*Ts*(zs+1))/...
 (4*m*(zs-1)^2+c2*(1+m/mdm)*Ts^2*(zs+1)^2+d2*(1+m/mdm)*2*(zs-1)*Ts*(zs+1));
opts.rho_min=0.85;
Delta_max_rho = ss(Delta_max); 
Delta_max_rho.A = Delta_max_rho.A/opts.rho_min;
Delta_max_rho.B = Delta_max_rho.B/opts.rho_min;
Dmax = hinfnorm(Delta_max_rho); % worst-case Hinfty norm
nu_iqc=2;
rhofilter=0.75;
P0 = [Dmax 0; 0 -1/Dmax];
iqc = IQC.dynamic(nu_iqc,rhofilter,P0);
%% controller synthesis
% parameters for optimization
opts.rho_start=0.99;
opts.iterations=5;
opts.rho_iter=10;
% performance channel specification
CH = struct();
CH.gam = -1; CH.gain = "p2p"; CH.z = 1:length(z); CH.w = 1:length(w);
CH.sigma = 0;
time_K = toc;
[Klarge,ana_ws,info] = iqc_synthesis(G_w,CH,iqc,opts);
gam = ana_ws.gam;
for i=1:length(info)
    time_K = time_K+info(i).time;
end
tic
% reduce order
K_reduce = reducespec(Klarge,"balanced");
K = getrom(K_reduce,Order = 2);
nK=length(K.A);
% Order = 2 yields no performance degradation as can be verified 
% by checking ana_ws.gam
info_ana = iqc_analysis(iqc,G_w,K,CH,opts);
gam_red = ana_ws.gam;
%% Observer synthesis
% set opts.only_K=true if no observer should be used and nu_t=1
% set opts.only_K=false to design an observer and optimize the initial 
% condition via nu_t
opts.only_K=false; % <<<<<<<<<<<<<<<<<<<<< use initial optimization?
opts.rho_start=ana_ws.rho;
Go = G_w([q z y y], :);
GK = lft(Go,K);
if ~opts.only_K
    [L,gamo1,~,~] = obs_synthesis(iqc,GK,info_ana.P,info_ana.M,info_ana.rho);
    % improve L by iterating between analysis and synthesis as we cannot 
    % simultaneously optimize L, P and M in a single SDP.
    [ana_obs] = obs_K_analysis(iqc,GK,L,opts);
    [L,~,~,~] = obs_synthesis(iqc,GK,ana_obs.P,ana_obs.M1+ana_obs.M2,ana_obs.rho);
    [ana_obs] = obs_K_analysis(iqc,GK,L,opts);
    gam2 = ana_obs.gam;
    [L,gamo2,muo,info] = obs_synthesis(iqc,GK,ana_obs.P,ana_obs.M1+ana_obs.M2,ana_obs.rho);
    time_L=toc;
else 
    L=ss([],[],[],zeros(nx+nK,ny));
end
%% Observer and Controller analysis
tic
param = struct();
param.rho=zeros(nz,1);
param.gam=zeros(nz,1);
param.mu=zeros(nz,1);
param.gamo=zeros(nz,1);
param.muo=zeros(nz,1);
param.PSchur = cell(nz,1);
param.cp0 = zeros(3,1);
param.cf0 = zeros(3,1);
npsi = iqc.n1+iqc.n2;
opts.rho_start = info.rho;
% Analyze each constraint direction individually as in Remark 2
for i=1:nz
    GKi = GK;
    GKi.OutputGroup.z = GK.OutputGroup.z(i);
    [ana_obs] = obs_K_analysis(iqc,GKi,L,opts);
    if ana_obs.info == "no solution found"
        error("no solution found")
    end
    param.rho(i) = ana_obs.rho;
    param.gam(i) = ana_obs.gam;
    param.mu(i) = ana_obs.mu;
    param.gamo(i) = ana_obs.gamo;
    param.muo(i) = ana_obs.muo;
    PXi11 = (ana_obs.P(1:npsi,1:npsi)-ana_obs.X1-ana_obs.X2);
    Pi12 = ana_obs.P(1:npsi,npsi+1:end);
    Pi22 = ana_obs.P(npsi+1:end,npsi+1:end);
    param.PSchur{i} = Pi22+Pi12'*(PXi11\Pi12);
    e1 = [0.005; -0.02];
    e2 = [0.005; 0.02];
    param.cf0(i) = e1'*Pi22(1:2,1:2)*e1;
    param.cf0(i) = max(param.cf0(i), e2'*Pi22(1:2,1:2)*e2);
    if ~opts.only_K
        Poi22 = ana_obs.Po(npsi+1:npsi+2,npsi+1:npsi+2);
        param.cp0(i) = e1'*Poi22*e1;
        param.cp0(i) = max(param.cp0(i), e2'*Poi22*e2);
    end
end
time_K_L_analysis = toc;
%% MPC simulation
tic
K.InputGroup.y=1:ny;
param.zmean = zmean;
param.zmin = zmin;
param.zmax = zmax;
param.dmax = 0.05/dscaling;
param.Q = 0.1*[K.D K.C 1]'*[K.D K.C 1]+blkdiag(100,1,1e-5*eye(length(K.A)+nu));
param.Q_to_design_T = eye(nx+length(K.A));
param.N = 45;
param.K = K;
param.L = L;
param.rf = 2;
param.only_K = opts.only_K;
opts_sim = struct();
opts_sim.Tsim = 60;
xe0 = [0.95; 0];
x0=xe0;
param.cp0=0*param.cp0;
param.cf0=0*param.cf0;
dsim = -param.dmax*ones(1,opts_sim.Tsim);
ADGc = [ 0                1               0     0     
         -c2*(1/m+1/mdm)  -d2*(1/m+1/mdm) c/m   d/m   
         0                0               0     1      
         c2/m             d2/m            -c/m  -d/m ];
BDGc = [0 0; 0 1/m; 0 0; 1/m 1/m];
CDG = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 1 0; 0 0 0 1];
DDG = [ 0 0; 0 0; 0 1; 0 0; 0 0 ];
ABDG = expm([ADGc BDGc; zeros(2,6)]*Ts);
ADGd = ABDG(1:4,1:4); 
BDGd = ABDG(1:4,5:6);
dDG = 1;
zDG = 1:3;
BDGd(:,dDG) = BDGd(:,dDG)*dscaling;
DDG(:,dDG) = DDG(:,dDG)*dscaling;
CDG(zDG,:) = diag(1./zscaling)*CDG(zDG,:);
DDG(zDG,:) = diag(1./zscaling)*DDG(zDG,:);
DG = lft(Delta_max, G_w);
xDG0 = [0;0;x0];
time_setup = toc;
[sim,mpc_data] = mpc_simulation(G_w,DG,xDG0,xe0,dsim,param,opts_sim);
if param.only_K
    time_offline = time_K+time_K_L_analysis+time_setup+sim.time_setup;
else
    time_offline = time_K+time_L+time_K_L_analysis+...
                   time_setup+sim.time_setup;
end
time_online = sim.time_sim/opts_sim.Tsim;
cost=0;
for t=1:50
cost = cost + sim.xDG(3:4,t)'*diag([100 1])*sim.xDG(3:4,t) + sim.u(:,t)'*0.1*sim.u(:,t);
end

%% plot
if opts.only_K==true
    disp("%%%%%%% RESULTS \nu_t=1 %%%%%%%")
    disp("gamma="+num2str(gam)+" after K synthesis.")
    %disp("gamma="+num2str(gam_red)+"after order reduction of K.")
    disp("gamma_1="+num2str(param.gam(1))+" after K analysis.")
    disp("gamma_2="+num2str(param.gam(2))+" after K analysis.")
    disp("gamma_3="+num2str(param.gam(3))+" after K analysis.")
    color = 'r';
else
    disp("%%%%%%% RESULTS optimize \nu_t %%%%%%%")
    disp("gamma="+num2str(gam)+" after K synthesis.")
    %disp("gamma="+num2str(gam_red)+"after order reduction of K.")
    disp("gamma^o="+num2str(gamo1)+" after L synthesis.")
    %disp("gamma="+num2str(gam2)+" after iterating L synthesis and K analysis.")
    %disp("gamma^o="+num2str(gamo2)+" after iterating L synthesis and K analysis.")
    disp("gamma_1="+num2str(param.gam(1))+" after simultaneous K & L analysis.")
    disp("gamma_2="+num2str(param.gam(2))+" after simultaneous K & L analysis.")
    disp("gamma_3="+num2str(param.gam(3))+" after simultaneous K & L analysis.")
    disp("gamma^0_1="+num2str(param.gamo(1))+" after simultaneous K & L analysis.")
    disp("gamma^0_2="+num2str(param.gamo(2))+" after simultaneous K & L analysis.")
    disp("gamma^0_3="+num2str(param.gamo(3))+" after simultaneous K & L analysis.")
    color = 'g';
end
disp("Closed loop cost until time Tsim: "+num2str(cost))
disp("x_(1,t) at t=43: "+num2str(sim.xDG(3,43)))
disp("Offline comp. time: "+time_offline+"s")
disp("Average online comp. time: "+time_online*1000+"ms")
figure(5)
plot(sim.xDG(3,:),sim.xDG(4,:),[color 'x'])      % real state
hold on
%plot(sim.xGKnf(1,:),sim.xGKnf(2,:),'bx') % nominal state
%plot(sim.xGKe(1,:),sim.xGKe(2,:),'rx')   % estimate
alpha = param.rho.^2./(1-param.rho.^2);
for t = 1:opts_sim.Tsim
    zt = sim.znf(:,t)-zmean(:);
    qdmax2 = sim.qnf(:,t).^2+param.dmax^2;
    ct = sqrt(param.gam./alpha.*sim.cf(:,t) + ...
              param.gam.*(param.gam-param.mu)*qdmax2);
    if length(ct)==1
        ct = ones(length(zt),1)*ct;
    end
    ztp = [sim.xGKe(1:2,t); sim.u(:,t)]./zscaling-zmean;
    qdmax2 = sim.qnp(:,t).^2+param.dmax^2;
    ctp = sqrt(param.gam./alpha.*sim.cp(:,t) + ...
    param.gam.*(param.gam-param.mu)*qdmax2);
    if length(ctp)==1
        ctp = ones(length(ztp),1)*ctp;
    end
    if mod(t,2)
    figure(5)
    ct = ct.*zscaling;
    plot(mpc_data(t).xGKnf(1,1)+[ct(1) ct(1) -ct(1) -ct(1) ct(1) ct(1)],...
        mpc_data(t).xGKnf(2,1)+[ct(2) -ct(2) -ct(2) ct(2) ct(2) -ct(2)],color)
    end
end