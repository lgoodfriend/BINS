% test BINS for convergence using Taylor Green as the exact solution
clear all
close all

advance
calculus

Lin = 2*pi;
BC= ones(4,1)*sqrt(-1);
nu = 0.1; % .................... molecular viscosity
IC_choice = 3;   % ...... initial condition
ng = 1; % ................... number of guardcells/ghost cells on each side of each dimension

idx=0;
for Ne=[20 30 40 50 60 70 80 90 100]
	idx=idx+1;
	% numerical solution
	Ne
	dx = Lin/Ne
	dt = dx*(1/100)
	BINS([Ne Ne],[Lin Lin],dt,0.1,BC,IC_choice,nu,ng);
	load( ['./BINS_output',int2str(Ne(1)),'.mat'], 'u','v', 'p', 'nu', 'L', 'N', 'T');
	
	% exact solution
	h=[dx dx];
	uexact = u; vexact = v;
	for i=ng+1:ng+N(1)+1; for j=ng+1:ng+N(2)+1
		% determine (x,y) coordinate of (i,j) point: left face
		x = (i-ng-1-0.5)*h(1); y = (j-ng-1)*h(2);
		uexact(i,j) = sin(x)*cos(y);
		% determine (x,y) coordinate of (i,j) point: bottom face
		x = (i-ng-1)*h(1); y = (j-ng-1-0.5)*h(2);
		vexact(i,j) = -cos(x)*sin(y);
	endfor; endfor;
	F = exp(-2*nu*T);
	uexact = uexact*F; vexact = vexact*F;
	
	% comparison
                        error(idx) = max(max(   abs(uexact(2:end-1,2:end-1)-u(2:end-1,2:end-1))    ));
                        deltaX(idx) = L(1)/N(1);                     
                        
                        clear XX YY x y uexact vexact F u v p L N T
endfor 

close all
% estimate slope
p = polyfit(log(deltaX),log(error),1);
order_of_convergence = p(1)
fit = exp(p(2))*deltaX.^p(1);

hold on
loglog(deltaX,error,'o')
loglog(deltaX,fit,'-')
hold off