function out = energyConstrainedPath(height, weight, Emax, Fmax, tfin, xnot, vnot, target)

%  Model parameters
params.m = weight; params.h = height; 
params.gamma = 0.2; % turning penalty

rho = 1.225;
Cd  = 0.9; 
Af = 0.0537985 * params.h^0.725 * params.m^0.425;
params.k = 0.5 * rho * Af * Cd;

params.aMax = 10.0; % not used
params.vEps = 0.10; % avoid division by zero

% TUNING param 
params.Emax = Emax;  % max allowed energy (z(T) <= Emax)
energyMargin = 0.10;         % 10% slack
params.Emin  = (1 - energyMargin) * Emax;
%  Target and timing
targetPos = target;
params.pTarget = targetPos; 
Tfinal  = tfin;
%Tfinal = 2.6;   

%  ICs
% px0 = 70.36;  py0 = 34.61;
% px1 = 69.93;  py1 = 34.21;
% vx0 = (px1 - px0)/ 0.1; vy0 = (py1 - py0)/ 0.1;
%vi_2 = [-2.76, -2.9];
% px0 = 72.83;  py0 = 36.33;
% vx0 = -1.54;  vy0 = -4.0;
px0 = xnot(1);  py0 = xnot(2);
vx0 = vnot(1);  vy0 = vnot(2);
%vx0 = -2.76;  vy0 =  -2.9;
z0  = 0.0;

x0  = [px0; py0; vx0; vy0; z0];

% OptimTraj problem setup
problem.func.dynamics = @(t,x,u) runnerDynamics_simple4(t,x,u,params);
problem.func.pathCst  = @(t,x,u) runnerPathCst2(t,x,u,params);
% minimize distance to target at final time
problem.func.bndObj = @(t0,x0,tF,xF) runnerTerminalCost2(t0,x0,tF,xF,targetPos);

% time bounds
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = Tfinal;
problem.bounds.finalTime.upp = Tfinal;

% state bounds: [px; py; vx; vy; z]
bounds.xLow  = [0; 0; -12; -12; 0]; 
bounds.xUpp  = [120; 53.3;  12;  12; 1e3];

problem.bounds.state.low = bounds.xLow;
problem.bounds.state.upp = bounds.xUpp;

% Initial state fixed to x0:
problem.bounds.initialState.low = x0;
problem.bounds.initialState.upp = x0;

% final state bounds
problem.bounds.finalState.low = bounds.xLow;
problem.bounds.finalState.upp = bounds.xUpp;
problem.bounds.finalState.low(5) = params.Emin;        
problem.bounds.finalState.upp(5) = params.Emax; 

% control bounds: u = [F; lambda]
problem.bounds.control.low = [0; 0];
problem.bounds.control.upp = [Fmax; 1];
%  initial guess
problem.guess.time = [0, Tfinal];
problem.guess.state = [x0, x0]; % start/end guess same

% initial guess for controls
F_init      = params.Emax / Tfinal;
lambda_init = 0.5;

problem.guess.control = [ ...
    F_init, F_init;   
    lambda_init, lambda_init];

%  solver options
problem.options.method = 'trapezoid';   % or 'trapezoid' 'hermiteSimpson' 'rungekutta4' 
problem.options.defaultAccuracy = 'low'; % very slow with medium
% manually set grid size
%problem.options.trapezoid.nGrid = tfin*10;
problem.options.nlpSolver = 'fmincon';

% problem.options.nlpOpt = optimoptions('fmincon', ...
%     'Algorithm', 'sqp', ... % try sqp
%     'Display', 'iter', ...
%     'MaxIterations', 5000, ...            
%     'MaxFunctionEvaluations', 5e5, ...
%     'OptimalityTolerance', 1e-4, ...       
%     'ConstraintTolerance', 1e-6, ...
%     'StepTolerance', 1e-8); 
problem.options.nlpOpt = optimset( ...
    'Display', 'iter-detailed', ...
    'MaxIter', 2000, ...
    'MaxFunEvals', 5e5, ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-6);

%  solve
soln = optimTraj(problem);

tSol = soln.grid.time;
xSol = soln.grid.state;
uSol = soln.grid.control;

pxSol = xSol(1,:); 
pySol = xSol(2,:); 
vxSol = xSol(3,:);
vySol = xSol(4,:); 
zSol = xSol(5,:); % accumulated energy

F_sol = uSol(1,:); % effort level
lambda_sol = uSol(2,:); % tradeoff

speedSol = sqrt(vxSol.^2 + vySol.^2);
out = struct('X', pxSol, 'Y', pySol, 'tSol', tSol, 'zSol', zSol, 'FSol', F_sol, 'lambdaSol', lambda_sol);
end