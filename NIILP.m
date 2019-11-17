%% Linear Program Solver by Y. Cui, K. Morikuni, T. Tsuchiya and K. Hayami.
% First version: Aug 2015
% Latest update: Nov 2019
% Please find the full article:
% Implementation of interior-point methods for LP based on 
% Krylov subspace iterative solvers with inner-iteration preconditioning
% Cui, Y., Morikuni, K., Tsuchiya, T. et al. 
% Comput Optim Appl (2019) 74: 143. 
% https://doi.org/10.1007/s10589-019-00103-y

%% main function: linear program interior point solver
function[x, z, info, trerrorHist] = NIILP(probname, solvername, tolmid, toladj, tol_lp, maxiter)
% The main function to solve linear programs.

% Args:
    % probname (string): the filename of the problem, e.g. 'lp_agg.mat'
    % tolmid (double): a tolerance level for the linear solver for the step equations
    % solvername (string): to specify linear solver, e.g. 'ABGMRES_scale'
    % toladj (bool): set to True to adjust the tolmid according to the progress
    % tol_lp (double): a tolerance level for the linear program solver
    % maxiter (int): max iteration number for the linear program solver
    
% Returns:
    % x (array of doubles): the primal solution
    % s (array of doubles): the dual solution
    % info (array): more information of this run, including:
    %     converged (double): 0 or 1 to indicate whether the linear program solver converged
    %     iter (double): number of iterations used for convergence
    %     avg1 (double): average number of iterations for the linear solver to solve the
    %           predictor step equations
    %     avg2 (double): average number of iterations for the linear solver to solve the
    %           corrector step equations
    % trerrorHist (array of doubles) the error at each interior-point step


% ------------------------------------------
% set up global variables
% ------------------------------------------
fprintf('\n<<<<< NIILP solver: %s for problem %s >>>>>\n', solvername, strtok(probname,'.'));
global x s y z w n m ub e Gt
global nub nt phi0 tau0 bnrm cnrm unrm
global Lbounds_non0 Ubounds_exist
global fixed_value_idx fixed_value

% ------------------------------------------
% setup interior-point parameters
% ------------------------------------------
tolmid_start = tolmid;
tau0 = .9995;       % step percentage to boundary
phi0 = 1.e-5;       % initial centrality factor

% ------------------------------------------
% load data & preprocess
% ------------------------------------------
[A, b, c, lb, ub] = prep(probname);
nt = n + nub;
e = ones(n,1);

% =============='''' tick ''''==============
t1 = cputime;

% initial point -----
bnrm =sqrt(b'*b); 
cnrm = sqrt(c'*c); 
unrm = [];
[x,y,z,s,w,Gt] = initp(A,b,c,ub,bnrm,cnrm); % Gt := A'y

if any(fixed_value_idx)
    fixed_value_idx(:);
    fixed_value(:);
end


if (Ubounds_exist) unrm = norm(ub); end

% ------------------------------------------
% compute the complementary quantities
% ------------------------------------------
Rxz = x.*z; Rsw = [];
if (Ubounds_exist) Rsw = s.*w; end
dgap = full(sum([Rxz; Rsw]));

% ------------------------------------------
% compute the feasibility measurements
% ------------------------------------------
[Rb, Rc, Ru, rb, rc, ru, rrb, rrc, rru, rdgap] = feas(A,b,c,dgap);
trerror = max([rrb rrc rru rdgap]);

% ------------------------------------------
% create empty variables to track progress
% ------------------------------------------
iter = 0; converged = 0; mzeros = sparse(m,1); nzeros = sparse(n,1);

fprintf('\n  Residuals:   Primal    Dual     Gap      TR_error\n');
fprintf('  ---------------------------------------------------------\n');

% ------------------------------------------
% interior point loops begin
% ------------------------------------------
while iter <= maxiter
    
    fprintf('  Iter %4i:  ', iter);
    fprintf('%8.2e %8.2e %8.2e %8.2e\n',rrb,rrc,rdgap,trerror);
    
    % ------------------------------------------
    % stop if convergence/divergence detected
    % ------------------------------------------
    if (iter > 0) && (trerror < tol_lp) && (trerror>=0)
        stop = 1; converged = 1; 
        if stop 
            fprintf('\n Successfully converged.\n')
            break
        end
    elseif iter > 0 && trerror < 0
        stop = 1; converged = 0;
        if stop 
            fprintf('\n Fail to converge: negative iterate ncountered.\n')
            break
        end
    end
    
    % ------------------------------------------
    % calculate the coefficient matrix for 
    % normal equations
    % ------------------------------------------
    xn1 = reciprocal(x);
    if (Ubounds_exist)
        sn1 = reciprocal(s);
        vmid = reciprocal(z.*xn1 + w.*sn1);
    else
        sn1 = []; 
        vmid = reciprocal(z.*xn1);
    end
    cap = 1.e+11; 
    vmid = full(min(cap,vmid));
    vmid_sqr = sqrt(vmid);
    G = 1./vmid_sqr;
    M = A*sparse(1:n,1:n,vmid_sqr,n,n,n);
    
    
    % ------------------------------------------
    % predictor stage
    % ------------------------------------------
    [dx,dz,ds,dw,solmid1(:,iter+1), Gdt, resmid1(iter+1), itmid1(iter+1), rhs1(:,iter+1)] = ...
        direct(A, M, Rb, Rc, Ru, Rxz, Rsw, vmid, G, xn1,sn1,z,w,0, tolmid, solvername);
    [ap,ad] = ratiotest(dx,dz,ds,dw);
    
    % ------------------------------------------
    % corrector stage
    % ------------------------------------------
    if (tau0*ap < 1 || tau0*ad < 1) % if we do need to take a corrector step
        mu = centering(dx,dz,ds,dw,ap,ad,dgap,trerror);
        Rxz = dx.*dz; Rsw = ds.*dw;
        [dx2,dz2,ds2,dw2,solmid2(:,iter+1), Gdt2, resmid2(iter+1), itmid2(iter+1), rhs2(:,iter+1)] = ...
            direct(A,M,mzeros,nzeros,nzeros, Rxz,Rsw,vmid,G, xn1,sn1,z,w,mu,tolmid,solvername);     
        dx = dx + dx2; 
        dz = dz + dz2;
        ds = ds + ds2;
        dw = dw + dw2;
        Gdt = Gdt+Gdt2; % Gdt := A'dy (saving 1 matrix-vector production)
        [ap,ad] = ratiotest(dx,dz,ds,dw);
    else % if the corrector step is not taken, we still need to put in some value for the record
        resmid2(iter+1) = 0;
        solmid2(:,iter+1) = zeros(n,1);
        itmid2(iter+1) = 0;
    end
    
    % ------------------------------------------
    % update iterate
    % ------------------------------------------
    [Rxz,Rsw,dgap,phi] = update(ap,ad,dx,dz,ds,dw, Gdt, trerror);  % phi: centrality measure
    iter = iter + 1;
    
    % ------------------------------------------
    % measure feasibility
    % ------------------------------------------
    [Rb, Rc, Ru, rb, rc, ru, rrb, rrc, rru, rdgap] = feas(A,b,c,dgap);
    trerror = max([rrb rrc rru rdgap]);
    
    % ------------------------------------------
    % adjust tolmid (optional)
    % ------------------------------------------
    if strcmp(toladj,'ON')
        if (resmid1(iter) > tolmid || resmid2(iter) > tolmid) && (tolmid < 1e-6)
            tolmid = 1.5*tolmid; % fail to converge
        end
        if log10(trerror)<= 1
            tolmid = max(1e-14,0.75*tolmid);
        end
        if log10(trerror)<= -3
            tolmid = max(1e-14, 0.5*tolmid);
        end
    end
    
    cn(iter) = cond(M*M');
    trerrorHist(iter) = trerror;
    
end
% ------------------------------------------
% interior point loops end
% ------------------------------------------

% =============='''' tick ends''''==============
t1 = cputime - t1;

% ------------------------------------------
% organise output
% ------------------------------------------
avg1 = round(sum(itmid1)/iter); avg2 = round(sum(itmid2)/iter);
info = [ converged; iter; avg1; avg2];

if (Lbounds_non0)
    x = x + lb; % put the lower bound back in
end

objp = full(c'*x);
info = [info; objp; t1];
fprintf('\n %s find the opt value: %8.2e \n', solvername, objp);
fprintf('\n  Iter  CPUtime  tolmid_start  tolmid_end   tolmidAdj \n');
fprintf('\n  %4i   %.2f     %8.2e         %8.2e          \n',...
             iter, t1,      tolmid_start, tolmid);

end

%% load and preprocess problem
function [A, b, c, lb, ub] = prep(probname)
% This function loads data and preprocess the problem.

% Args:
    % probname (string): the filename of the problem, e.g. 'lp_agg.mat'
    
% Returns:
    % the matrices and arrays of the problem

global Lbounds_non0
global Ubounds_exist nub m n
global fixed_value_idx fixed_value

% ------------------------------------------
% load problem
% ------------------------------------------
problem_loaded = load(probname);
A = problem_loaded.Problem.A; 
b = problem_loaded.Problem.b; 
c = problem_loaded.Problem.aux.c; 
lb = problem_loaded.Problem.aux.lo; 
ub = problem_loaded.Problem.aux.hi;

if any(lb > ub)
    fprintf('\n Error in NIILP.prep(): problem %s has lower bound exceeding upper bound\n', probname); return
end


% ------------------------------------------
% delete zero rows in A
% ------------------------------------------
indices = sum(A~=0, 2)>0; 
b = b(indices);
A = A(indices, :);

fixed_value_idx = find(lb==ub & lb == 0);
fixed_value = ub(fixed_value_idx);

[m,n] = size(A);


% ------------------------------------------
% shift nonzero lower bounds
% ------------------------------------------
Lbounds_non0 = any(lb ~= 0);
if (Lbounds_non0)
    lb(isinf(lb))= -1e-11;
    b = b - A*lb; % when lb has -inf, need to adjust b to be a big number
end


% ------------------------------------------
% find noninf upper bounds
% ------------------------------------------
nub = 0;
position = isfinite(ub); % mark the finite bounds with 1
Ubounds_exist = any(position);

ub(~isfinite(ub))=0; 
ub = position.*ub; 
if (Ubounds_exist)
    nub = nnz(ub);
    ub = sparse(position.*(ub-lb));
end

end

%% choose an initial guess
function [x, y, z, s, w, Gt] = initp(A, b, c, ub, bnrm, cnrm)
% This function sets up a initial guess for the interior-point process
% following the LIPSOL and the Lustig 1992 paper.

global Ubounds_exist e
global fixed_value_idx fixed_value


[m,n] = size(A); y = zeros(m,1);
pmin = max(bnrm/100, 100);
dmin = cnrm*.425; dmin = max(dmin, pmin/40);
pmin = min(pmin, 1.e+4); dmin = min(dmin, 1.e+3);
P = A*sparse(1:n,1:n,e,n,n,n)*A';
rho = min(100,bnrm);

solspd = P\(full(b-rho*A*e)); 

x = A'*solspd + rho*e;
pmin = max(pmin, -min(x)); 
x = max(pmin,x);

x(fixed_value_idx) = fixed_value;


z = full((c+dmin).*(c > 0) + dmin*(-dmin < c & c <= 0) - c.*(c <= -dmin));
s = []; w = [];

if (Ubounds_exist)
    s = spones(ub).*max(pmin, ub-x);
    w = spones(ub).*( dmin*(c > 0) + ...
        (dmin - c).*(-dmin < c & c <= 0) - 2*c.*(c <= -dmin) );
end
Gt = A'*y;

end

%% calculate feasibility
function [Rb, Rc, Ru, rb, rc, ru, rrb, rrc, rru, rdgap] = feas(A, b, c, dgap)
% This function returns the (relative) feasibility measure and duality measure.

global x z s w Gt
global ub bnrm cnrm unrm nt
global Ubounds_exist


% ------------------------------------------
% compute vectors
% ------------------------------------------
Rb = A*x - b;
Rc = Gt + z - c;
Ru = [];
if (Ubounds_exist)
    Rc = Rc - w;
    Ru = spones(s).*x + s - ub;
end

% ------------------------------------------
% compute norms
% ------------------------------------------
rb = norm(Rb); rc = norm(Rc); ru = 0;
if Ubounds_exist 
    ru = norm(Ru); 
end

if any(isnan([rb rc ru])) 
    error(' NaN occured.'); 
end

% ------------------------------------------
% compute relative error
% ------------------------------------------
rrb = rb/max(1,bnrm); rrc = rc/max(1,cnrm); rru = 0;
if (Ubounds_exist) 
    rru = ru/(1+unrm);  
end
rdgap = dgap/nt;

end

%% calculate direction to move
function [dx, dz, ds, dw, dt, Gdt, resmid, itmid, rhs] = ...
direct(A, M, Rb, Rc, Ru, Rxz, Rsw, vmid, G, xn1, sn1, z, w, mu, tolmid, solvername)
% This function computes the direction by calling a pre-specified linear algebera solver. 


% ------------------------------------------
% construct rhs
% ------------------------------------------
global Ubounds_exist n m
maxitmid = m;

if (mu ~= 0) 
    Rxz = Rxz - mu; 
end

Rc = Rc - Rxz.*xn1;

if (Ubounds_exist)
    if (mu ~= 0) 
        Rsw = Rsw - mu; 
    end
    Rc = Rc + (Rsw - Ru.*w).*sn1;
end

rhs = -(Rb + A*(vmid.*Rc));


% ------------------------------------------
% select a solver
% ------------------------------------------
if strcmp(solvername,'ABGMRES_scale')
    [dt, reshis, itmid] = ABNESOR4IP_scale(M', rhs, tolmid);
    resmid = reshis(itmid);
    Gdt = G.* dt;        
elseif strcmp(solvername,'CGNE_scale')
    [dt, reshis, itmid] = CGNE4IP_scale(M', rhs, tolmid, maxitmid);
    resmid = reshis(itmid);
    Gdt = G.* dt;  
elseif strcmp(solvername,'MRNE_scale')
    [dt, reshis, itmid] = MRNE4IP_scale(M', rhs, tolmid, maxit);
    resmid = reshis(itmid);
    Gdt = G.* dt; 
elseif strcmp(solvername,'modDirect')
    [dy,resmid] = moDirect(M, rhs);
    itmid = 0; dt = zeros(n,1);
    Gdt = A'*dy;
else
    fprintf('Please check your input: undefined solver')
end


dx = vmid.*(Gdt + Rc);
dz = -(z.*dx + Rxz).*xn1;
ds = []; dw = [];

if (Ubounds_exist)
    ds = -(dx.*spones(w) + Ru);
    dw = -(w.*ds + Rsw).*sn1;
end

end

%% choose centering parameter
function mu = centering(dx,dz,ds,dw,ap,ad,dgap,trerror)
% This function calculates a centering parameter.

global x z s w
global Ubounds_exist nt

newdgap = (x + min(1,ap)*dx)'*(z + min(1,ad)*dz);
if (Ubounds_exist)
    newdgap = newdgap + (s + min(1,ap)*ds)'*(w + min(1,ad)*dw);
end
sigmak = (newdgap/dgap)^2;
sigmin = 0; sigmax = .208;

p = ceil(log10(trerror));
if (p < -2 && dgap < 1.e+3) sigmax = 10^(p+1); end

sigmak = max(sigmin, min(sigmax, sigmak));

mu = sigmak*dgap/nt;

end

%% choose a step length
function [ap,ad] = ratiotest(dx,dz,ds,dw)
% This function finds a step length.
global x z s w
global Ubounds_exist

% ------------------------------------------
% ratio test
% ------------------------------------------
ap = -1/min([dx./x; -0.1]);
ad = -1/min([dz./z; -0.1]);
if (Ubounds_exist)
    as = -1/min([ds(find(s))./nonzeros(s); -0.1]);
    aw = -1/min([dw(find(w))./nonzeros(w); -0.1]);
    ap = min(ap, as); ad = min(ad, aw);
end

end

%% update the iterate
function [Rxz,Rsw,dgap,phi] = update(ap,ad,dx,dz,ds,dw, Gdt, trerror)
% This function update the iterate and 
% calculate the duality measure and centrality measure.

global x z s w Gt n ub
global phi0 tau0 % phi0 is set at 1e-5; tau0 is set at 0.9995;
global backed nt
global Ubounds_exist

tau = tau0;
if ~exist('backed') backed = 0; end
if ~backed tau = .9 + 0.1*tau0; end
k = ceil(log10(trerror));
if (k <= -5) tau = max(tau,1-10^k); end

ap = min(tau*ap,1); ad = min(tau*ad,1);
xc = x; %yc = y;
Gtc = Gt;
zc = z;
sc = s;
wc = w;

step = [1 .9975 .95 .90 .75 .50];
for k = 1:length(step)
    x = xc + step(k)*ap*dx;
    Gt = Gtc + step(k)*ad*Gdt;
    z = zc + step(k)*ad*dz;
    s = sc + step(k)*ap*ds;
    w = wc + step(k)*ad*dw;
    Rxz = x.*z; Rsw = [];
    if (Ubounds_exist) Rsw = s.*w; end
    
    dgap = full(sum([Rxz; Rsw]));
    
    swtmp = Rsw(ub~=0);
    phi = nt*full(min([Rxz; swtmp]))/dgap;

    
    if (max(ap,ad) == 1) || (phi >= phi0)  
        break; 
    end

end

phi0 = min(phi0, phi); if (k > 1) backed = 1; end

end

%% helper function: reciprocal with cap
function Y = reciprocal(X)
% This is a helper function to ceil the reciprocal.

if issparse(X)
    [m, n]  = size(X);
    [i,j,Y] = find(X);
    Y = sparse(i,j,1./Y,m,n);
else
    Y = 1./X;
end

ceiling = 1.e+16; 

Y = min(ceiling,Y);

end
% --- end of this file ---
