% test accuracy of BINS by comparing to exact solution to Taylor-Green vortex
clear all
close all

advance
calculus

L = 2*pi; % length of domain
N = 70; % number of points in domain
BC= ones(4,1)*sqrt(-1); % boundary conditions
nu = 0.1; % molecular viscosity
IC_choice = 3;   % initial condition
ng = 1; % number of guardcells/ghost cells on each side of each dimension
T = 0.1; % end time


% numerical solution
dx = L/N
dt = dx*(1/10) % grid size determines time step for stability
BINS(N,L,dt,T,BC,IC_choice,nu,ng);
load( ['./BINS_output',int2str(N),'.mat'], 'u','v', 'p', 'nu', 'L', 'N', 'T');

% exact solution
h=dx;
uexact = u; vexact = v;
for i=ng+1:ng+N+1; for j=ng+1:ng+N+1
	% determine (x,y) coordinate of (i,j) point: left face
	x = (i-ng-1)*h; y = (j-ng-1 + 0.5)*h;
	uexact(i,j) = sin(x)*cos(y);
	% determine (x,y) coordinate of (i,j) point: bottom face
	x = (i-ng-1 + 0.5)*h; y = (j-ng-1)*h;
	vexact(i,j) = -cos(x)*sin(y);
endfor; endfor;
F = exp(-2*nu*T);
uexact = uexact*F; vexact = vexact*F;

% compare numerical to exact solution
u_error = abs( u-uexact );
v_error = abs( v-vexact );

plotSoln(u_error,v_error,ng,L,N,h)