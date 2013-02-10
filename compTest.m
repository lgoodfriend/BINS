% test accuracy of BINS by comparing to exact solution to Taylor-Green vortex
clear all
close all

advance
calculus

% choose solution parameters--------------------------------------------------------------------------
L = 2*pi; % length of domain
N = 40; % number of points in domain
BC= ones(4,1)*sqrt(-1); % boundary conditions
nu = 0.1; % molecular viscosity
IC_choice = 3;   % initial condition
ng = 1; % number of ghost cells on each side of each dimension

dx = L/N;
dt = dx*(1/10); % grid size determines time step for stability
T = 1; % end time

% numerical solution---------------------------------------------------------------------------------------
BINS(N,L,dt,T,BC,IC_choice,nu,ng);
load( ['./BINS_output',int2str(N),'.mat'], 'u','v', 'p', 'nu', 'L', 'N', 'T');

% exact solution---------------------------------------------------------------------------------------------
h=dx;
uexact = u; vexact = v;
li = ng+1; % lower index of inner (non-ghost cell) data
ui = ng+N; % upper index of inner (non-ghost cell) data
% for u
%
xx = linspace(0,L,N+1);
yy = linspace(h/2,L-h/2,N);
[y,x] = meshgrid(yy,xx);
uexact(li:ui+1,li:ui) = sin(x).*cos(y);
% for v
%
xx = linspace(h/2,L-h/2,N);
yy = linspace(0,L,N+1);
[y,x] = meshgrid(yy,xx);
vexact(li:ui,li:ui+1) = -cos(x).*sin(y);
% decay with time
F = exp(-2*nu*T);
uexact = uexact*F; vexact = vexact*F;

% compare numerical to exact solution and plot-----------------------------------------------
u_error = abs( u-uexact );
v_error = abs( v-vexact );

plotSoln(u_error,v_error,ng,L,N,h)