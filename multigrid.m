% DANGER
% VERY BUGGY
% DO NOT USE


% functions for solving Poisson pressure equation Ap = f
% using multigrid with red-black ordered Gauss-Seidel
% contains the functions
%	multigrid_solve
%	mg
%	laplacian
%	gauss_seidel
%	average
%	interpolate
1; % this is needed in octave for the interpreter to find all the functions
%------------------------------------------------------------------------------------------------------------------------------------------
function p = multigrid_solve(f,h,N,BC)
	% solves laplacian(p) = f
	% inputs:
	% f: right hand side of Poisson equation (N x N array)
	% h: grid spacing dx and dy
	% N: number of grid points in x and y
	% BC: array defining boundary conditions (4 x 1 array)
	% 
	% returns:
	% p: the solution to the Poisson equation, the new pressure
	p = zeros(size(f));
	error = 100000;
	comp = 1;
	iter = 0;
	not_done = true;
	while not_done % loop through iterations of the multigrid cycle
		% one iteration
		p = mg(p,f,h,N,BC);
		iter = iter+1;
		% are we done?
		LaP = laplacian(p,h,N,BC);
		error = max(max(abs( LaP -f )));
		comp =0.001*max(max(abs( f )));
		% if we're done, return
		iter
		error
		comp
		if (error<comp); return; endif;
		if (iter>100); return; endif;
	endwhile
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
function p = mg(p,f,h,N,BC)
	% recursive function that goes up and down the multigrid "V"
	% inputs:
	% p: the current solution to the Poisson equation (N x N array)
	% f: the current right hand side to the Poisson equation (N x N array)
	% h: the current grid spacing dx and dy
	% N: the current number of grid points in x and y
	% BC: the array of boundary conditions (4 x 1 array)
	% 
	% returns:
	% p: updated solution to the Poisson equation
	
	% relax solution with Gauss-Seidel
	p = gauss_seidel(p,f,h,N,BC);
	p = gauss_seidel(p,f,h,N,BC);
	if (N>1);
		% calculate remainder
		LaP = laplacian(p,h,N,BC);
		r = LaP - f;
		% go down the "V"
		rc = average(r,N);
		delta = zeros(size(rc));
		rc = mg(delta,rc,2*h,N/2,BC);
		% go up the "V"
		r = interpolate(delta,N);
		p = p+r;
	endif
	p = gauss_seidel(p,f,h,N,BC);
	p = gauss_seidel(p,f,h,N,BC);	
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
function LaP = laplacian(p,h,N,BC)
	% make domain with ghost cells
	pGC = zeros(N+2,N+2);
	pGC(2:N+1,2:N+1) = p;
	% fill domain's ghost cells
	foo = zeros(N+3,N+3); % needed as u and v placeholders
	[foo,foo,pGC] = fillBC(foo,foo,pGC,1,N,BC);
	
	% calculate laplacian on inner domain
	LaP = zeros(size(p));
	LaP = (pGC(1:N,2:N+1) - 2*pGC(2:N+1,2:N+1) + pGC(3:N+2,2:N+1))/(h^2) ...% d^2p/dx^2
	+ (pGC(2:N+1,1:N) - 2*pGC(2:N+1,2:N+1) + pGC(2:N+1,3:N+2))/(h^2); % d^2p/dy^2
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
function p = gauss_seidel(p,f,h,N,BC)
	% make domain with ghost cells
	pGC = zeros(N+2,N+2);
	pGC(2:N+1,2:N+1) = p;
	% fill domain's ghost cells
	foo = zeros(N+3,N+3); % needed as u and v placeholders
	[foo,foo,pGC] = fillBC(foo,foo,pGC,1,N,BC);
	
	% relax solution
	% red squares
	pGC(2:2:N,2:2:N) = 0.25*( pGC(1:2:N-1,2:2:N) + pGC(3:2:N+1,2:2:N) ...% left + right
		+ pGC(2:2:N,1:2:N-1) + pGC(2:2:N,3:2:N+1) ...% + down + up
		- (h^2)*f(1:2:N-1,1:2:N-1) ); % right hand side
	pGC(3:2:N+1,3:2:N+1) = 0.25*( pGC(2:2:N,3:2:N+1) + pGC(4:2:N+2,3:2:N+1) ...% left + right
		+ pGC(3:2:N+1,2:2:N) + pGC(3:2:N+1,4:2:N+2) ...% + down + up
		- (h^2)*f(2:2:N,2:2:N) ); % right hand side
	p(1:2:N-1,1:2:N-1) = pGC(2:2:N,2:2:N);
	p(2:2:N,2:2:N) = pGC(3:2:N+1,3:2:N+1);
	[foo,foo,pGC] = fillBC(foo,foo,pGC,1,N,BC);
	% black squares
	p(2:2:N,1:2:N-1) = 0.25*( pGC(2:2:N,2:2:N) + pGC(4:2:N+2,2:2:N) ...% left + right
		+ pGC(3:2:N+1,1:2:N-1) + pGC(3:2:N+1,3:2:N+1) ...% + down + up
		- (h^2)*f(2:2:N,1:2:N-1) ); % right hand side
	p(1:2:N-1,2:2:N) = 0.25*( pGC(1:2:N-1,3:2:N+1) + pGC(3:2:N+1,3:2:N+1) ...% left + right
		+ pGC(2:2:N,2:2:N) + pGC(2:2:N,4:2:N+2) ...% + down + up
		- (h^2)*f(1:2:N-1,2:2:N) ); % right hand side
	
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
function rc = average(r,N)
	% average solutions in a 2x2 region to approximate the solution at the next
	% smaller grid
	rc = zeros(N/2,N/2);
	rc = 0.25*( r(2:2:N,2:2:N) + r(1:2:N-1,2:2:N) + r(2:2:N,1:2:N-1) + r(1:2:N-1,1:2:N-1) );
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
function r = interpolate(delta,N)
	% inject solutions from the coarse to the fine grid
	r = zeros(N,N);
	r(2:2:N,2:2:N) = delta(1:N/2,1:N/2);
	r(1:2:N-1,2:2:N) = delta(1:N/2,1:N/2);
	r(2:2:N,1:2:N-1) = delta(1:N/2,1:N/2);
	r(1:2:N-1,1:2:N-1) = delta(1:N/2,1:N/2);
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------








