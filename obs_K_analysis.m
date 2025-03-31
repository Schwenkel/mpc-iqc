function [ana_obs] = obs_K_analysis(iqc,GK,L,options)
if nargin == 3
    options = struct();
end
if ~isfield(options,"rho_min")    % lower bound on rho for which the IQC
    options.rho_min = 0;          % is still valid for Delta_rho
end
if ~isfield(options,"margin_pos") % margin for strict inequalities, i.e.,
    options.margin_pos = 5e-8;    % LMI > 0   <==>   LMI >= margin_pos
end
if ~isfield(options,"only_K")     % margin for strict inequalities, i.e.,
    options.only_K = false;       % LMI > 0   <==>   LMI >= margin_pos
end
if ~isfield(options,"rho_start")  % start rho optimization 
    options.rho_start = 0.99;
end
options.rho_iter = 15;    % LMI > 0   <==>   LMI >= margin_pos
eps = options.margin_pos;
sdp_opts = sdpsettings('solver','mosek','verbose',0);
% channels and dimensions:
rho_iter=options.rho_iter;
rho = options.rho_start;
[LMIs,margins,ana_var] = get_LMIs(iqc,GK,L,rho,eps, options.only_K);
w_gam = 2;
optimize(LMIs, w_gam*ana_var.gam+ana_var.gamo, sdp_opts);
ana_obs = get_values(ana_var);
error_occured = false;
if check_LMIs(LMIs, margins) < 0
    error_occured=true;
    gam = Inf;
else
    gam = w_gam*ana_obs.gam+ana_obs.gamo;
end
persistent rhotry
if isempty(rhotry)
    rhotry = rho-0.01;
end
for i=1:rho_iter
    [LMIs,margins,ana_var] = get_LMIs(iqc,GK,L,rhotry,eps,options.only_K);
    assign_warm_start(ana_var, ana_obs);
    optimize(LMIs, w_gam*ana_var.gam+ana_var.gamo, sdp_opts);
    if check_LMIs(LMIs, margins)<0
        gamtry = Inf;
    else
        error_occured = false;
        ana_ws_try = get_values(ana_var);
        gamtry = w_gam*ana_ws_try.gam+ana_ws_try.gamo;
        if gamtry < gam  % best gama so far
            ana_obs = ana_ws_try;
        end
    end
    [rhotry, rho, gam] = step_size_rule(rhotry,gamtry,rho,gam,options.rho_min);
end
if error_occured
    ana_obs.info = "no solution found";
else
    ana_obs.info = "solution found";
end
end

% Auxiliary functions:

function [LMIs, margins, ana_var] = get_LMIs(iqc, GK, L, rho, margin, only_K)
    if only_K
        [LMIs, margins, ana_var] = get_LMIs_only_K(iqc, GK, rho, margin);
        return
    end
    % loop trafo:
    GKrho = GK;
    GKrho.A = GKrho.A/rho;
    GKrho.B = GKrho.B/rho;
    % Building the system Sigma:
    S = build_Sigma(iqc.Psi1,iqc.Psi2,GKrho);
    % Dimensions and channels
    nP = size(S.A,1);
    % SDP variables
    P_ = sdpvar(nP);
    gam_ = sdpvar(1);
    mu_ = sdpvar(1);
    ana_var = struct();
    ana_var.P = P_;
    ana_var.mu = mu_;  % mu = gam whenever not p2p
    ana_var.H = iqc.H;
    ana_var.Psi1 = iqc.Psi1;
    ana_var.Psi2 = iqc.Psi2;
    ana_var.rho = rho;
    Lrho = L;
    Lrho.A = L.A/rho;
    Lrho.B = L.B/rho;
    Xi = build_Xi(S, Lrho);
    % Dimensions and channels
    nPo = size(Xi.A,1);
    Po_  = sdpvar(nPo);
    gamo_ = sdpvar(1);
    muo_  = sdpvar(1);
    M1 = iqc.M;
    X1 = iqc.X;
    [M2, X2, MX2_con, H2] = iqc.get_MX2;
    M = M1+M2;
    X = X1+X2;
    [M3, X3, MX3_con, H3] = iqc.get_MX2;
    [M4, X4, MX4_con, H4] = iqc.get_MX2;
    H = [iqc.H; H2; H3; H4];
    % SDP constraints
    LMIs = [iqc.MX_con{1}, MX2_con{1}, MX3_con{1}, MX4_con{1}, gamo_>=muo_, muo_>=0];
    LMIs = [LMIs LMI1(Xi, M+M3+M4, muo_, Po_, margin)];
    LMIs = [LMIs LMI2(Xi, M4,X3,X4,muo_,gamo_,Po_,P_,rho,margin)];
    margins = [iqc.MX_con{2}; MX2_con{2};MX3_con{2};MX4_con{2}; 0; 0; margin; margin];
    ana_var.Po = Po_;
    ana_var.M1 = M1;
    ana_var.M2 = M2;
    ana_var.X1 = X1;
    ana_var.X2 = X2;
    ana_var.M3 = M3;
    ana_var.M4 = M4;
    ana_var.X3 = X3;
    ana_var.X4 = X4;
    ana_var.gamo = gamo_;
    ana_var.muo = muo_;
    ana_var.H = H;
    alpha = rho^2/(1-rho^2);
    beta = mu_;
    margins = [margins; margin; 0.01*margin];
    p2pLMI1 = pe2p_LMI1(S, P_, M, mu_, margins(end-1));
    p2pLMI2 = pe2p_LMI2(alpha, beta, S, P_, X1, X2, M2, gam_, margins(end));
    LMIs = [LMIs, p2pLMI1, p2pLMI2];
    margins = [margins; 10*margin];
    LMIs = [LMIs, PX_pos_LMI(P_, X, margins(end))];
    ana_var.gam = gam_;
end

function [LMIs, margins, ana_var] = get_LMIs_only_K(iqc, GK, rho, margin)
    % loop trafo:
    GKrho = GK;
    GKrho.A = GKrho.A/rho;
    GKrho.B = GKrho.B/rho;
    % Building the system Sigma:
    S = build_Sigma(iqc.Psi1,iqc.Psi2,GKrho);
    % Dimensions and channels
    nP = size(S.A,1);
    % SDP variables
    P_ = sdpvar(nP);
    gam_ = sdpvar(1);
    mu_ = sdpvar(1);
    ana_var = struct();
    ana_var.P = P_;
    ana_var.mu = mu_;  % mu = gam whenever not p2p
    ana_var.H = iqc.H;
    ana_var.Psi1 = iqc.Psi1;
    ana_var.Psi2 = iqc.Psi2;
    ana_var.rho = rho;
    % Dimensions and channels
    M1 = iqc.M;
    X1 = iqc.X;
    [M2, X2, MX2_con, H2] = iqc.get_MX2;
    M = M1+M2;
    X = X1+X2;
    H = [iqc.H; H2];
    % SDP constraints
    LMIs = [iqc.MX_con{1}, MX2_con{1}];
    margins = [iqc.MX_con{2}; MX2_con{2}];
    ana_var.Po = 1;
    ana_var.M1 = M1;
    ana_var.M2 = M2;
    ana_var.X1 = X1;
    ana_var.X2 = X2;
    ana_var.M3 = 0;
    ana_var.M4 = 0;
    ana_var.X3 = 0;
    ana_var.X4 = 0;
    ana_var.gamo = 0;
    ana_var.muo = 0;
    ana_var.H = H;
    alpha = rho^2/(1-rho^2);
    beta = mu_;
    margins = [margins; margin; 0.01*margin];
    p2pLMI1 = pe2p_LMI1(S, P_, M, mu_, margins(end-1));
    p2pLMI2 = pe2p_LMI2(alpha, beta, S, P_, X1, X2, M2, gam_, margins(end));
    LMIs = [LMIs, p2pLMI1, p2pLMI2];
    margins = [margins; 10*margin];
    LMIs = [LMIs, PX_pos_LMI(P_, X, margins(end))];
    ana_var.gam = gam_;
end

function LMI = PX_pos_LMI(P, X, margin)
    nP = length(P);
    LMI = P - blkdiag(X,zeros(nP-length(X))) >= margin*eye(nP);
end

function LMI = pe2p_LMI1(S, P, M, mu, margin)
    ns = length(M);
    np = length(S.InputGroup.p);
    s = 1:ns;
    nwi = size(S.B,2)-np;
    nP = size(P,1);
    outer = [eye(nP)          zeros(nP, np+nwi);
             S.A              S.B
             S.C(s,:)         S.D(s,:)
             zeros(nwi, nP)   zeros(nwi,np)  eye(nwi)];
    inner = blkdiag(-P, P, M, -mu*eye(nwi));
    LM = outer'*inner*outer;
    LMI = LM <= -margin*eye(length(LM));
end

function LMI = pe2p_LMI2(alp, bet, S, P, X1, X2, M2, gam, margin)
    ns = length(M2);
    nP = length(P);
    s = 1:ns;
    np = length(S.InputGroup.p);
    z = S.OutputGroup.z;  nzi = length(z);
    nwi = size(S.B,2)-np;
    X1bar = blkdiag(X1,zeros(nP-length(X1)));
    X2bar = blkdiag(X2,zeros(nP-length(X1)));
    if isfloat(X2) ...            % if D_Sigma^zp == 0 and sigma=0
        && isfloat(M2) ...        % then we get rid of the p row/col as it 
        && min(min(X2==0)) ...    % is all zeros 
        && min(min(M2==0)) ...
        && min(min(S.D(z,p)==0)) % then:
        inner = blkdiag(X1bar-P, -alp*(gam-bet)*eye(nw));
        % Schur complement
        w = S.InputGroup.w;
        LM = [inner    [S.C(z,:)'; S.D(z,w)'];
              S.C(z,:)  S.D(z,w)    -gam/alp*eye(nz)];
    else
        outer = [eye(nP)          zeros(nP,np+nwi)
                 S.A              S.B
                 S.C(s,:)         S.D(s,:)
                 zeros(nwi,nP+np)  eye(nwi)];
        inner = blkdiag(X1bar-P, X2bar, M2, -alp*(gam-bet)*eye(nwi));
        % Schur complement
        LM = [outer'*inner*outer      [S.C(z,:)'; S.D(z,:)'];
              S.C(z,:)  S.D(z,:)    -gam/alp*eye(nzi)];
    end
    LM = LM + LM'; % numeric noise can cause LM to be non-symmetric
    LMI = LM <= -2*margin*eye(length(LM));
end

function assign_warm_start(ana_var, ana_ws)
    if isfield(ana_ws, "P")
        if length(ana_ws.P) == length(ana_var.P)
            assign(ana_var.P,   ana_ws.P)
        end
        assign(ana_var.mu,  ana_ws.mu)
        assign(ana_var.gam, ana_ws.gam)
        if isfield (ana_var, "M") 
            if ~isfloat(ana_var.X) % If float and no sdpvar "assign" produces an error
                assign(ana_var.X, ana_ws.X)
            end
            if ~isfloat(ana_var.M) % If float and no sdpvar "assign" produces an error
                assign(ana_var.M, ana_ws.M)
            end
            if ~isempty(ana_var.H)
                assign(ana_var.H(:), ana_ws.H(1:length(ana_var.H)));
            end
        else
            if ~isfloat(ana_var.X1) % If float and no sdpvar "assign" produces an error
                assign(ana_var.X1, ana_ws.X1)
                assign(ana_var.X2, ana_ws.X2)
            end
            if ~isfloat(ana_var.M1) % If float and no sdpvar "assign" produces an error
                assign(ana_var.M1, ana_ws.M1)
                assign(ana_var.M2, ana_ws.M2)
            end
            if ~isempty(ana_var.H)
                if length(ana_var.H) == length(ana_ws.H)
                    assign(ana_var.H, ana_ws.H);
                else
                    assign(ana_var.H(:), [ana_ws.H(:); ana_ws.H(:)]);
                end
            end
        end
    end
end

function ana_ws = get_values(ana_var)
    ana_ws = struct();
    if ~isempty(ana_var.H)
        ana_ws.H = value(ana_var.H);
    else
        ana_ws.H = [];
    end
    ana_ws.P    = value(ana_var.P);
    ana_ws.M1   = value(ana_var.M1);
    ana_ws.X1   = value(ana_var.X1);
    ana_ws.M2   = value(ana_var.M2);
    ana_ws.X2   = value(ana_var.X2);
    ana_ws.mu   = value(ana_var.mu);
    ana_ws.gam  = value(ana_var.gam);
    ana_ws.Psi1 = ana_var.Psi1;
    ana_ws.Psi2 = ana_var.Psi2;
    ana_ws.rho  = ana_var.rho;
    ana_ws.Po   = value(ana_var.Po);
    ana_ws.M3   = value(ana_var.M3);
    ana_ws.M4   = value(ana_var.M4);
    ana_ws.X3   = value(ana_var.X3);
    ana_ws.X4   = value(ana_var.X4);
    ana_ws.gamo = value(ana_var.gamo);
    ana_ws.muo  = value(ana_var.muo);
    ana_ws.H    = value(ana_var.H);
end

function feasibility = check_LMIs(LMIs, margins)
    % if feasibility > 0 then feasible
    tol = 1e-7 ;
    feasibility = min(check(LMIs)+margins*4/5+tol);
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

function Sigma = build_Sigma(Psi1,Psi2,GKrho)
    p = GKrho.InputGroup.p;
    w = GKrho.InputGroup.w;
    q = GKrho.OutputGroup.q;
    z = GKrho.OutputGroup.z;
    y = GKrho.OutputGroup.y;
    n1 = length(Psi1.A);
    n2 = length(Psi2.A);
    nx = length(GKrho.A);
    nw = length(w);
    np = length(p);
    ny = length(y);
    nz = length(z);
    A = [Psi1.A          zeros(n1,n2)  Psi1.B*GKrho.C(q,:);
         zeros(n2,n1)    Psi2.A        zeros(n2,nx)
         zeros(nx,n1)    zeros(nx,n2)  GKrho.A];
    B = [ Psi1.B*GKrho.D(q,p)  Psi1.B*GKrho.D(q,w); 
          Psi2.B               zeros(n2,nw);
          GKrho.B(:,p)         GKrho.B(:,w)];
    C = [ Psi1.C               Psi2.C         Psi1.D*GKrho.C(q,:)
          zeros(nz,n1)         zeros(nz,n2)   GKrho.C(z,:)
          zeros(ny,n1)         zeros(ny,n2)   GKrho.C(y,:)];
    D = [ Psi2.D+Psi1.D*GKrho.D(q,p)  Psi1.D*GKrho.D(q,w)
          GKrho.D(z,p)                GKrho.D(z,w)
          GKrho.D(y,p)                GKrho.D(y,w)];
    Sigma = ss(A,B,C,D,-1);
    Sigma.InputGroup.p = 1:np;
    Sigma.InputGroup.w = np+1:np+nw;
    ns = size(Psi1.C,1);
    Sigma.OutputGroup.s = 1:ns;
    Sigma.OutputGroup.z = ns+1:ns+nz;
    Sigma.OutputGroup.y = ns+nz+1:ns+nz+ny;
end

function Xi = build_Xi(Sigma,Lrho)
    p = Sigma.InputGroup.p;
    w = Sigma.InputGroup.w;
    s = Sigma.OutputGroup.s;
    y = Sigma.OutputGroup.y;
    nL = length(Lrho.A);
    nPhi = length(Sigma.A);
    nx = size(Lrho.C,1);
    ns = length(s);
    nw = length(w);
    np = length(p);
    A = [ Sigma.A              zeros(nL); 
          Lrho.B*Sigma.C(y,:)  Lrho.A ];
    B = [ Sigma.B; Lrho.B*Sigma.D(y,:)];
    C = [ Sigma.C(s,:)    zeros(ns,nL)
          eye(nPhi-nx)  zeros(nPhi-nx,nx+nL)
          [zeros(nx,nPhi-nx) eye(nx)]-Lrho.D*Sigma.C(y,:)  -Lrho.C];
    D = [ Sigma.D(s,:); zeros(nPhi-nx,np+nw); -Lrho.D*Sigma.D(y,:) ];
    Xi = ss(A,B,C,D,-1);
    Xi.InputGroup.p = 1:np;
    Xi.InputGroup.w = np+1:np+nw;
    Xi.OutputGroup.s = 1:ns;
    Xi.OutputGroup.z = ns+1:ns+nL;
end


function [rhotry, rhonew, gam] = step_size_rule(rhotry,gamtry,rho,gam,rhomin)
    if gamtry<gam  % direction was good, go 2x in that direction
        rhonew = rhotry;
        gam = gamtry;
        rhotry = rhotry+2*(rhotry-rho);
    else % direction was bad, go 0.33x in opposite direction
        rhonew = rho;
        rhotry = rho+(rho-rhotry)/3;
    end
    if abs(rhotry-rhonew) < 1e-4 % make sure steps do not get too small
        rhotry = rhonew+sign(rhotry-rhonew)*5e-2;
    end
    if rhotry>=1 % ensure rho in (rhomin,1)
        rhotry = 1/2+rhonew/2;
    elseif rhotry<=rhomin
        rhotry = rhomin/2 + rhonew/2;
    end
end
