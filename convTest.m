% test BINS order of convergence using Taylor-Green as the exact solution
clear all
close all

% choose solution parameters--------------------------------------------------------------------------
Lin = 2*pi; % length of domain
BC = [sqrt(-1) sqrt(-1) sqrt(-1) sqrt(-1)]; % boundary conditions
nu = 0.1; % molecular viscosity
IC_choice = 3;   % initial condition: Taylor-Green vortex
ng = 1; % number of ghost cells on each side of each dimension
Tin = 2; % end time

idx=0;
% loop over a range of grid sizes-----------------------------------------------------------------------
for Ne=[40 50 60 70 80 90 100 110 120] 
	idx=idx+1;
	% numerical solution
	Ne
	dx = Lin/Ne
	dt = dx*(1/10) % grid size determines time step for stability
	BINS(Ne,Lin,dt,Tin,BC,IC_choice,nu,ng);
	load( ['./BINS_output',int2str(Ne),'.mat'], 'u','v', 'p', 'ng','nu', 'L', 'N', 'T');
	
	% exact solution
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
	% include decay with time
	%
	F = exp(-2*nu*T)
	uexact = uexact*F; vexact = vexact*F;
	
	% compare numerical to exact solution
	uerr =  max(max(   abs(uexact(2:end-1,2:end-1)-u(2:end-1,2:end-1))    ));
	verr =  max(max(   abs(vexact(2:end-1,2:end-1)-v(2:end-1,2:end-1))    ));
	error(idx) = max(  [uerr verr]  );
	deltaX(idx) = L/N;                     
                         
	clear XX YY x y uexact vexact F u v p L N T
end
close all

% calculate convergence rate and plot---------------------------------------------------------------
% estimate slope using linear fit
p = polyfit(log(deltaX),log(error),1);
order_of_convergence = p(1)
fit = exp(p(2))*deltaX.^p(1);

% plot!
loglog(deltaX,error,'o','MarkerSize',10)
hold on
loglog(deltaX,fit,'-','LineWidth',2)
hold off
xlabel('\Delta x')
ylabel('Error')
grid on
xlim([5/120 7/40])
ylim([.003 .02])

save( ['./conv_data.mat'], ...
	'deltaX','error','fit');
