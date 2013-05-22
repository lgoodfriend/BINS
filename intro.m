% introduce BINS and its use
% Taylor-Green vortex example
close all
clear all

% choose parameters for this solution
N = 32; % N = number of cells in x and y
L = 2*pi; % L = physical dimension in x and y
dt = (L/N)/10; % dt = timestep
T = 1; % T = final time
BC = [sqrt(-1) sqrt(-1) sqrt(-1) sqrt(-1)]; % BC = boundary conditions (4x1 array)
IC_choice = 3; % IC_choice = choice of initial conditions defined in IC.m: Taylor-Green vortex
nu = 0.1; % nu = molecular viscosity
ng = 1; % ng = number of ghostcells 

% solve NS
BINS(N,L,dt,T,BC,IC_choice,nu,ng)
