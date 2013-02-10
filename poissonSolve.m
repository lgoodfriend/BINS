function solution = poissonSolve(rhs,h,BC,ng,N)
	% solve the Poisson equation:   d^2 solution / dx^2 + d^2 solution/dy^2 = rhs
	% Numerically, this equation results in a matrix inversion:
	% A_ij solution_j = rhs_i
	% for more details, see Section 3.4 on page 8
	%
	% two solution methods are implemented here: 
	%	matlab's built-in matrix inverter
	%	multigrid with a red-black Gauss Seidel core (not done yet!)
	%
	% inputs:
	% rhs: right hand side to Poisson equation, size [N+2*ng   N+2*ng]
	% h: grid size
	% BC: array defining boundary conditions, size [4]
	% ng: number of guardcells
	% N: number of points in the grid
	%
	% returns:
	% solution: the solution to the Poisson equation, the new pressure, size [N+2*ng   N+2*ng]
	solution = rhs;
	
	% solve using multigrid (not done yet!)
	%solution(ng+1:ng+N,ng+1:ng+N) = multigrid_solve(ng+1:ng+N,ng+1:ng+N),h,N,BC);
	%return
	
	% solve using built in matrix inverter
	% make a lookup table for converting between matrix and vector representations
	% of the xy plane
	bigM = N^2;
	for i=1:N; for j=1:N
		m(i,j) = i + (j-1)*N;
		Mlookup( m(i,j), :) = [i,j];                    	          	
	endfor; endfor
          	                          
	% make the matrix A that defines the Poisson equation with the BCs
	% and the matrix f that defines our right hand side          	                    	
	[A,f] = make_matrix(bigM,N,ng,h,Mlookup,m,rhs,BC);
                                	           
	% solve solnVector = A\f using matlab's built-ins
	solnVector = A\f;
                                	           
	% convert back to x,y,z coordinates
	for Midx = 1:bigM
		idx = Mlookup(Midx,:);
		i = idx(1); j=idx(2);
		solution(i+ng,j+ng) = solnVector(Midx);
	endfor

endfunction







