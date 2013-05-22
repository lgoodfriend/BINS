% test accuracy of BINS by comparing to exact solution to Taylor-Green vortex
clear all
close all

% choose solution parameters------------------------------------------------------------------------
L = 2*pi; % length of domain
N = 64; % number of points in domain
BC= ones(4,1)*sqrt(-1); % boundary conditions
nu = 0.1; % molecular viscosity
IC_choice = 3;   % initial condition: Taylor-Green vortex
ng = 1; % number of ghost cells on each side of each dimension

dx = L/N;
dt = dx*(1/10); % grid size determines time step for stability
T = 2; % end time

% numerical solution--------------------------------------------------------------------------------
BINS(N,L,dt,T,BC,IC_choice,nu,ng);
load( ['./BINS_output',int2str(N),'.mat'], 'u','v', 'p', 'nu', 'L', 'N', 'T');

% exact solution------------------------------------------------------------------------------------
h=dx;
uexact = u; vexact = v;
li = ng+1; % lower index of inner (non-ghost cell) data
ui = ng+N; % upper index of inner (non-ghost cell) data
% for u
xu = linspace(0,L,N+1);
yu = linspace(h/2,L-h/2,N);
[y,x] = meshgrid(yu,xu);
uexact(li:ui+1,li:ui) = sin(x).*cos(y);
% for v
xv = linspace(h/2,L-h/2,N);
yv = linspace(0,L,N+1);
[y,x] = meshgrid(yv,xv);
vexact(li:ui,li:ui+1) = -cos(x).*sin(y);
% decay with time
F = exp(-2*nu*T);
uexact = uexact*F; vexact = vexact*F;

% compare numerical to exact solution and plot------------------------------------------------------
% pick out line y=pi/2 in exact and numerical solutions
y_line = pi/2;
yIdx = round(y_line/L * N); 
u_line = u(li:ui+1,yIdx); v_line = v(li:ui,yIdx);
ue_line = uexact(li:ui+1,yIdx); ve_line = vexact(li:ui,yIdx);
% plot
figure()
hold on
plot(xu,ue_line,'k-',xu,u_line,'ko')
plot(xv,ve_line,'b-',xv,v_line,'bo')
hold off
xlabel('x'); ylabel('velocity');
legend('u: exact','u: numerical','v: exact','v: numerical')
title_str = sprintf('Taylor Green solution at y = %f',y_line); title(title_str)
xlim([0 L]);
