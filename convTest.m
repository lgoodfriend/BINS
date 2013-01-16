% test BINS for convergence using Taylor Green as the exact solution
clear all
close all

advance
calculus

Lin = 2*pi; % length of domain
BC= ones(4,1)*sqrt(-1); % boundary conditions
nu = 0.1; % molecular viscosity
IC_choice = 3;   % initial condition
ng = 1; % number of guardcells/ghost cells on each side of each dimension
Tin = 0.1; % end time

idx=0;
% loop over a range of grid sizes
for Ne=[70 80 90 100] 
	idx=idx+1;
	% numerical solution
	Ne
	dx = Lin/Ne
	dt = dx*(1/10) % grid size determines time step for stability
	BINS(Ne,Lin,dt,Tin,BC,IC_choice,nu,ng);
	load( ['./BINS_output',int2str(Ne),'.mat'], 'u','v', 'p', 'nu', 'L', 'N', 'T');
	
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
                        error(idx) = max(max(   abs(uexact(2:end-1,2:end-1)-u(2:end-1,2:end-1))    ));
                        deltaX(idx) = L/N;                     
                        
                        clear XX YY x y uexact vexact F u v p L N T
endfor 

close all
% estimate slope
p = polyfit(log(deltaX),log(error),1);
order_of_convergence = p(1)
fit = exp(p(2))*deltaX.^p(1);

% plot!
hold on
loglog(deltaX,error,'o','MarkerSize',24)
loglog(deltaX,fit,'-','LineWidth',2)
hold off