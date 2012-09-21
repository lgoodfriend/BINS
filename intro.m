% introduce BINS and its use

% open access to functions defined in larger files
advance
calculus 

% choose parameters for this solution
N = [20 20]; % N = number of cells in x and y (2x1 array)
L = [2*pi 2*pi]; % L = physical dimension in x and y (2x1 array)
dt = 0.003; % dt = timestep
T = 10*dt; % T = final time
BC = [sqrt(-1) sqrt(-1) sqrt(-1) sqrt(-1)]; % BC = boundary conditions (4x1 array)
IC_choice = 3; % IC_choice = choice of initial conditions defined in IC.m
nu = 0.01; % nu = molecular viscosity
ng = 1; % ng = number of ghostcells 

% solve NS
BINS(N,L,dt,T,BC,IC_choice,nu,ng)