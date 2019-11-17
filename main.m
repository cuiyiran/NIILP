
%% Linear Program Solver by Y. Cui, K. Morikuni, T. Tsuchiya and K. Hayami.
% First version: Aug 2015 
% Latest update: Nov 2019
% Please find the full article:
% Implementation of interior-point methods for LP based on 
% Krylov subspace iterative solvers with inner-iteration preconditioning
% Cui, Y., Morikuni, K., Tsuchiya, T., Hayami, K.
% Comput Optim Appl (2019) 74: 143. 
% <https://doi.org/10.1007/s10589-019-00103-y>

%% Main Script
% This is the main script to run the linear programming solver.
% The functions are in a separate file NIILP.m.
% Please edit the configurable fields below to change the embedded linear solvername,
% tolerance, or problem to solve.

clear;
clc;
warning off

% ---------------------------------------------------
% choose a problem
% this is the file name of the problem
% other options: see more LPnetlib problems 
% in <https://sparse.tamu.edu>
% ---------------------------------------------------
probname = 'lp_adlittle.mat';

% ---------------------------------------------------
% choose a linear solver
% this is the file name of the solver
% other options:
% 'CGNE_scale'
% 'MRNE_scale'
% 'modDirect' (to be released)
% ---------------------------------------------------
solvername = 'ABGMRES_scale';

% ---------------------------------------------------
% choose a (starting) tolerance for the linear solvers 
% ---------------------------------------------------
tolmid = 1e-6;

% ---------------------------------------------------
% choose whether to adjust the linear solver tolerance 
% as the interior-point iterations progress
% ---------------------------------------------------
toladj='ON';

% ---------------------------------------------------
% choose a tolerance for the linear program solver
% ---------------------------------------------------
tol_lp = 1e-8;

% ---------------------------------------------------
% choose a max number of iterations for the linear 
% program solver
% ---------------------------------------------------
max_iter = 100;

% ---------------------------------------------------
% run the solver
% ---------------------------------------------------

[x, z, info, error_vec] = NIILP(probname, solvername, tolmid, toladj, tol_lp, max_iter);

% ---------------------------------------------------
% simple plot of the interior-point progress
% ---------------------------------------------------

if info(1)==1
    hold off
    hold on
    grid on
    plot(log10(error_vec(:)), '-*','MarkerSize',5);
    
    legend('interior point error');
    hx = xlabel('interior-point iteration step');
    hy = ylabel('log10 of the interior point error');
   
    title(sprintf('progress of NIILP solving problem %s with solver %s', strtok(probname,'.'), solvername),'Interpreter', 'none');

end
