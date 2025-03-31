function [L,gamo,muo,info] = obs_synthesis(iqc,GK,P,M,rho_start)
eps = 5e-8; % margin for strict inequalities
options = sdpsettings('solver','mosek','verbose',0);
% channels and dimensions:
rho_iter=0;
rho = rho_start;
[LMIs, margins, variables] = get_LMIs(iqc,GK,P,M,rho,eps);
optimize(LMIs, variables.gamo, options);
obs_values = get_values(variables);
error_occured = false;
if check_LMIs(LMIs, margins) < 0
    error_occured=true;
    agamo = Inf;
    rhotry = rho/4+3/4;
else
    agamo = rho^2/(1-rho^2)*obs_values.gamo;
    rhotry = rho-0.01;
end
for i=1:rho_iter
    [LMIs, margins, variables] = get_LMIs(iqc,GK,P,M,rhotry,eps);
    optimize(LMIs, variables.gamo, options);
    if check_LMIs(LMIs, margins)<0
        agamtry = Inf;
    else
        error_occured = false;
        obs_val_try = get_values(variables);
        agamtry = rhotry^2/(1-rhotry^2)*obs_val_try.gamo;
        if agamtry < agamo  % best gama so far
            obs_values = obs_val_try;
        end
    end
    [rhotry, rho, agamo] = step_size_rule(rhotry,agamtry,rho,agamo);
    if rhotry>=1
        [rhotry, rho, agamo] = step_size_rule(rhotry,inf,rho,agamo);
    end
end
L = obs_values.Lrho;
L.A = L.A*rho;
L.B = L.B*rho;
muo = obs_values.muo;
info = obs_values;
info.rho = rho;
gamo = info.gamo;
if error_occured
    info.info="No solution found. ";
else
    info.info = "solution found";
end
end

% Auxiliary functions:

function [LMIs, margins, variables] = get_LMIs(iqc,GK,P,M,rho,margin)
    % loop trafo:
    GKrho = GK;
    GKrho.A = GKrho.A/rho;
    GKrho.B = GKrho.B/rho;
    % Building the system Phi
    Phi = build_Phi(iqc.Psi1, iqc.Psi2, GKrho);
    % Dimensions and channels
    nPo = size(Phi.A,1);
    nx = length(GKrho.A);
    y = Phi.OutputGroup.y;
    ny = length(y);
    P1o_  = sdpvar(nPo);
    P2o_  = sdpvar(nPo);
    gamo_ = sdpvar(1);
    muo_  = sdpvar(1);
    Ko_   = sdpvar(nPo);
    Lo_   = sdpvar(nPo, ny);
    Mo_   = sdpvar(nx,nPo);
    No_   = sdpvar(nx,ny);
    ABo_  = [ Lo_*Phi.C(y,:)+Ko_  Lo_*Phi.C(y,:)  Lo_*Phi.D(y,:)];
    J2 = [zeros(nx,iqc.n1+iqc.n2) eye(nx)];
    CDo_  = [ J2-No_*Phi.C(y,:)-Mo_  J2-No_*Phi.C(y,:)  -No_*Phi.D(y,:)];
    % SDP constraints
    M3 = iqc.M;
    X3 = iqc.X;
    [M4, X4, MX4_con, H4] = iqc.get_MX2;
    H = [iqc.H; H4];
    LMIs = [iqc.MX_con{1}, MX4_con{1}, gamo_>=muo_, muo_>=0, ...
            P1o_>=margin*eye(nPo), P2o_>=margin*eye(nPo)];
    margins = [iqc.MX_con{2}; MX4_con{2}; 0; 0; margin; margin; margin; margin];
    LMIs = [LMIs LMI1(Phi, ABo_, M+M3+M4, muo_, P1o_, P2o_, margin)];
    LMIs = [LMIs LMI2(Phi,CDo_,M4,X3,X4,muo_,gamo_,P1o_,P2o_,P,rho,margin)];
    variables=struct();
    variables.P1o = P1o_;
    variables.P2o = P2o_;
    variables.M3 = M3;
    variables.M4 = M4;
    variables.X3 = X3;
    variables.X4 = X4;
    variables.Ko = Ko_;
    variables.Lo = Lo_;
    variables.Mo = Mo_;
    variables.No = No_;
    variables.gamo = gamo_;
    variables.muo = muo_;
    variables.H = H;
end

function LMI = LMI1(Phi, ABo_, Mo_, muo_, P1o_, P2o_, margin)
    nPo = length(P1o_);
    s = Phi.OutputGroup.s;
    np = length(Phi.InputGroup.p);
    nw = length(Phi.InputGroup.w);
    LM = -blkdiag([P1o_ P1o_; P1o_ P2o_], zeros(np+nw));
    AB = [ Phi.A  Phi.A  Phi.B];
    ABABo_ = AB'*ABo_;
    LM = LM + ABABo_ + ABABo_' + AB'*P2o_*AB;
    CD = [ Phi.C(s,:) Phi.C(s,:) Phi.D(s,:)];
    LM = LM + CD'*Mo_*CD;
    LM = LM - muo_*blkdiag(zeros(2*nPo+np), eye(nw));
    LM = [LM ABo_'; ABo_ P1o_-P2o_];
    LMI = LM+LM' <= -2*margin*eye(length(LM));
end

function LMI = LMI2(Phi,CDo_,M4_,X3_,X4_,muo_,gamo_,P1o_,P2o_,P,rho,margin)
    nPo = length(P1o_);
    npsi = length(X4_);
    nx = length(Phi.A)-npsi;
    s = Phi.OutputGroup.s;
    np = length(Phi.InputGroup.p);
    nw = length(Phi.InputGroup.w);
    X4_ = blkdiag(X4_, zeros(nPo-npsi));
    X3_ = blkdiag(X3_, zeros(nPo-npsi));
    LM = -blkdiag([P1o_-X3_ P1o_-X3_; P1o_-X3_ P2o_], zeros(np+nw));
    AB = [ Phi.A  Phi.A  Phi.B];
    LM = LM + AB'*X4_*AB;
    CD = [ Phi.C(s,:) Phi.C(s,:) Phi.D(s,:)];
    LM = LM + CD'*M4_*CD;
    J1 = [eye(npsi) zeros(npsi,nx)];
    P11 = P(1:npsi,1:npsi);
    P12 = P(1:npsi,npsi+1:end);
    P22 = P(npsi+1:end,npsi+1:end);
    CDo1 = [J1 J1 zeros(npsi,nw+np)];
    CDoPCDo_ = CDo1'*P12*CDo_;
    LM = LM + CDoPCDo_' + CDoPCDo_ + CDo1'*P11*CDo1;
    LM = LM - rho^2/(1-rho^2)*(gamo_-muo_)*blkdiag(zeros(2*nPo+np), eye(nw));
    LM = [LM CDo_'; CDo_ -inv(P22)];
    LMI = LM <= -margin*eye(length(LM));
end

function obs_values = get_values(variables)
    nPo = length(variables.P1o);
    Po_diff = value(variables.P2o-variables.P1o);
    Po = [value(variables.P2o) eye(nPo); eye(nPo) inv(Po_diff)];
    AL = -value(variables.Ko)/Po_diff;
    BL = value(variables.Lo);
    CL = -value(variables.Mo)/Po_diff;
    DL = value(variables.No);
    obs_values=struct();
    obs_values.Lrho = ss(AL,BL,CL,DL,-1);
    obs_values.Po   = Po;
    obs_values.M3   = value(variables.M3);
    obs_values.M4   = value(variables.M4);
    obs_values.X3   = value(variables.X3);
    obs_values.X4   = value(variables.X4);
    obs_values.gamo = value(variables.gamo);
    obs_values.muo  = value(variables.muo);
    obs_values.H    = value(variables.H);
end

function feasibility = check_LMIs(LMIs, margins)
    % if feasibility > 0 then feasible
    tol = 1e-7 ;
    feasibility = min(check(LMIs)+margins*4/5+tol);
end

function Phi = build_Phi(Psi1,Psi2,GKrho)
    p = GKrho.InputGroup.p;
    w = GKrho.InputGroup.w;
    q = GKrho.OutputGroup.q;
    y = GKrho.OutputGroup.y;
    n1 = length(Psi1.A);
    n2 = length(Psi2.A);
    nx = length(GKrho.A);
    nw = length(w);
    np = length(p);
    ny = length(y);
    A = [Psi1.A zeros(n1,n2) Psi1.B*GKrho.C(q,:);
         zeros(n2,n1) Psi2.A zeros(n2,nx)
         zeros(nx,n1+n2) GKrho.A];
    B = [ Psi1.B*GKrho.D(q,p)  Psi1.B*GKrho.D(q,w); 
          Psi2.B               zeros(n2,nw);
          GKrho.B(:,p)         GKrho.B(:,w)];
    C = [ Psi1.C               Psi2.C         Psi1.D*GKrho.C(q,:)
          zeros(ny,n1)         zeros(ny,n2)   GKrho.C(y,:)];
    D = [ Psi2.D+Psi1.D*GKrho.D(q,p)  Psi1.D*GKrho.D(q,w)
          GKrho.D(y,p)         GKrho.D(y,w)];
    Phi = ss(A,B,C,D,-1);
    Phi.InputGroup.p = 1:np;
    Phi.InputGroup.w = np+1:np+nw;
    ns = size(Psi1.C,1);
    Phi.OutputGroup.s = 1:ns;
    Phi.OutputGroup.y = ns+1:ns+ny;
end
