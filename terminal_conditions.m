function [S,T,cff,GKxnff] = terminal_conditions(Q,AGK,CGKq,CGKz,param)
nz = length(param.rho);
nGK = length(AGK);
S = dlyap(AGK',Q(1:nGK,1:nGK));
T = chol(dlyap(AGK',param.Q_to_design_T));
cff = param.mu.*param.rho.^2./(1-param.rho.^2).*param.dmax^2;
param.rf = min(1,param.rf); % param.rf must be >= 1 such that (18c) holds
GKxnff = zeros(nz,1); 
a = sqrt(max(sum(CGKz/T.^2))); 
for i=1:nz
    b = param.gam(i)^2*param.dmax^2;
    c = param.gam(i)*(param.gam(i)+(param.rf-1)*param.mu(i))*norm(CGKq/T)^2;
    % maximize GKxnff s.t. a*sqrt(GKxnff)+sqrt(b+c*GKxnff) <= 1 via
    % bisection
    GKxnff_min = 0;
    GKxnff_max = 1/(a+sqrt(c))^2;
    for j=1:10
        GKxnff_try = (GKxnff_min+GKxnff_max)/2;
        if a*sqrt(GKxnff_try)+sqrt(b+c*GKxnff_try) <= 1
            GKxnff_min = GKxnff_try;
        else
            GKxnff_max = GKxnff_try;
        end
    end
    GKxnff(i) = GKxnff_min;
    cff(i) = param.mu(i)*param.rho(i)^2/(1-param.rho(i)^2)*...
                      (param.dmax^2+param.rf*norm(CGKq/T)^2*GKxnff(i));
end
GKxnff = min(GKxnff); % terminal constraint is the smallest.
% For the other constraints, cff may be increased:
cff = (1-a*sqrt(GKxnff))^2-param.gam.*(param.gam-param.mu)*norm(CGKq/T)^2*GKxnff;
cff = param.rho.^2./(1-param.rho.^2)./param.gam.*cff;
end



