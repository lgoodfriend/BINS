% introduce BINS and its use
% cavity flow example
close all
clear all

% choose parameters for this solution
N = 32; % N = number of cells in x and y
L = 1; % L = physical dimension in x and y
dt = (L/N) / 10; % dt = timestep
T = 0.1;%10*dt; % T = final time
BC = [0 0 0 1]; % BC = boundary conditions
IC_choice = 1; % IC_choice = choice of initial conditions defined in IC.m
nu = 0.01; % nu = molecular viscosity
ng = 1; % ng = number of ghostcells 

% solve NS
BINS(N,L,dt,T,BC,IC_choice,nu,ng)
