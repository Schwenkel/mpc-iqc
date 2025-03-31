function [gamo,muo,info] = obs_analysis(iqc,GK,L,ana_ws_all)
eps = 5e-8; % margin for strict inequalities
options = sdpsettings('solver','mosek','verbose',0);
% channels and dimensions:
for i=1:length(ana_ws_all)
    P=ana_ws_all{i}.P;
    M=ana_ws_all{i}.M;
    rho=ana_ws_all{i}.rho;
    [LMIs, margins, variables] = get_LMIs(iqc,GK,L,P,M,rho,eps);
    sol = optimize(LMIs, variables.gamo, options);
    obs_values = get_values(variables);
    if check_LMIs(LMIs, margins) < 0
        info(i) = "No solution found. ";
    else
        info(i) = "solution found";
    end
    muo(i) = obs_values.muo;
    gamo(i) = obs_values.gamo;
end
end

% Auxiliary functions:

function [LMIs, margins, variables] = get_LMIs(iqc,GK,L,P,M,rho,margin)
    % loop trafo:
    GKrho = GK;
    GKrho.A = GKrho.A/rho;
    GKrho.B = GKrho.B/rho;
    % Building the system Phi
    Phi = build_Phi(iqc.Psi1, iqc.Psi2, GKrho);
    Lrho = L;
    Lrho.A = L.A/rho;
    Lrho.B = L.B/rho;
    Xi = build_Xi(Phi, Lrho);
    % Dimensions and channels
    nPo = size(Xi.A,1);
    Po_  = sdpvar(nPo);
    gamo_ = sdpvar(1);
    muo_  = sdpvar(1);
    M3 = iqc.M;
    X3 = iqc.X;
    [M4, X4, MX4_con, H4] = iqc.get_MX2;
    H = [iqc.H; H4];
    LMIs = [iqc.MX_con{1}, MX4_con{1}, gamo_>=muo_, muo_>=0];
    LMIs = [LMIs LMI1(Xi, M+M3+M4, muo_, Po_, margin)];
    LMIs = [LMIs LMI2(Xi, M4,X3,X4,muo_,gamo_,Po_,P,rho,margin)];
    margins = [iqc.MX_con{2}; MX4_con{2}; 0; 0; margin; margin];
    variables=struct();
    variables.Po = Po_;
    variables.M3 = M3;
    variables.M4 = M4;
    variables.X3 = X3;
    variables.X4 = X4;
    variables.gamo = gamo_;
    variables.muo = muo_;
    variables.H = H;
end

function LMI = LMI1(Xi, Mo_, muo_, Po_, margin)
    nPo = length(Po_);
    s = Xi.OutputGroup.s;
    np = length(Xi.InputGroup.p);
    nw = length(Xi.InputGroup.w);
    inner = blkdiag(-Po_, Po_, Mo_, -muo_*eye(nw));
    outer = [ eye(nPo)       zeros(nPo,np+nw);
              Xi.A           Xi.B; 
              Xi.C(s,:)      Xi.D(s,:)
              zeros(nw,nPo)  zeros(nw,np) eye(nw)];
    LM = outer'*inner*outer;
    LM = (LM+LM')/2;
    LMI = LM <= -margin*eye(length(LM));
end

function LMI = LMI2(Xi,M4_,X3_,X4_,muo_,gamo_,Po_,P,rho,margin)
    nPo = length(Po_);
    npsi = length(X4_);
    np = length(Xi.InputGroup.p);
    nw = length(Xi.InputGroup.w);
    X4_ = blkdiag(X4_, zeros(nPo-npsi));
    X3_ = blkdiag(X3_, zeros(nPo-npsi));
    beta = rho^2/(1-rho^2)*(gamo_-muo_);
    inner = blkdiag(X3_-Po_, X4_, M4_, P, -beta*eye(nw));
    outer = [ eye(nPo)       zeros(nPo,np+nw);
              Xi.A           Xi.B; 
              Xi.C           Xi.D
              zeros(nw,nPo)  zeros(nw,np) eye(nw)];
    LM = outer'*inner*outer;
    LMI = LM <= -margin*eye(length(LM));
end

function obs_values = get_values(variables)
    obs_values=struct();
    obs_values.Po   = value(variables.Po);
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
    A = [Psi1.A          zeros(n1,n2)  Psi1.B*GKrho.C(q,:);
         zeros(n2,n1)    Psi2.A        zeros(n2,nx)
         zeros(nx,n1)    zeros(nx,n2)  GKrho.A];
    B = [ Psi1.B*GKrho.D(q,p)  Psi1.B*GKrho.D(q,w); 
          Psi2.B               zeros(n2,nw);
          GKrho.B(:,p)         GKrho.B(:,w)];
    C = [ Psi1.C               Psi2.C         Psi1.D*GKrho.C(q,:)
          zeros(ny,n1)         zeros(ny,n2)   GKrho.C(y,:)];
    D = [ Psi2.D+Psi1.D*GKrho.D(q,p)  Psi1.D*GKrho.D(q,w)
          GKrho.D(y,p)                GKrho.D(y,w)];
    Phi = ss(A,B,C,D,-1);
    Phi.InputGroup.p = 1:np;
    Phi.InputGroup.w = np+1:np+nw;
    ns = size(Psi1.C,1);
    Phi.OutputGroup.s = 1:ns;
    Phi.OutputGroup.y = ns+1:ns+ny;
end

function Xi = build_Xi(Phi,Lrho)
    p = Phi.InputGroup.p;
    w = Phi.InputGroup.w;
    s = Phi.OutputGroup.s;
    y = Phi.OutputGroup.y;
    nL = length(Lrho.A);
    nPhi = length(Phi.A);
    nx = size(Lrho.C,1);
    ns = length(s);
    nw = length(w);
    np = length(p);
    A = [ Phi.A              zeros(nL); 
          Lrho.B*Phi.C(y,:)  Lrho.A ];
    B = [ Phi.B; Lrho.B*Phi.D(y,:)];
    C = [ Phi.C(s,:)    zeros(ns,nL)
          eye(nPhi-nx)  zeros(nPhi-nx,nx+nL)
          [zeros(nx,nPhi-nx) eye(nx)]-Lrho.D*Phi.C(y,:)  -Lrho.C];
    D = [ Phi.D(s,:); zeros(nPhi-nx,np+nw); -Lrho.D*Phi.D(y,:) ];
    Xi = ss(A,B,C,D,-1);
    Xi.InputGroup.p = 1:np;
    Xi.InputGroup.w = np+1:np+nw;
    Xi.OutputGroup.s = 1:ns;
    Xi.OutputGroup.z = ns+1:ns+nL;
end