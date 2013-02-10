% introduce BINS and its use
% Taylor-Green vortex example

% open access to functions defined in larger files
advance
calculus 

% choose parameters for this solution
N = 32; % N = number of cells in x and y
L = 2*pi; % L = physical dimension in x and y
dt = (L/N)/10; % dt = timestep
T = 1; % T = final time
BC = [sqrt(-1) sqrt(-1) sqrt(-1) sqrt(-1)]; % BC = boundary conditions (4x1 array)
IC_choice = 3; % IC_choice = choice of initial conditions defined in IC.m
nu = 0.1; % nu = molecular viscosity
ng = 1; % ng = number of ghostcells 

% solve NS
BINS(N,L,dt,T,BC,IC_choice,nu,ng)