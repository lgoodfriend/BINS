function solution = poissonSolve(rhs,h,BC,ng,N,LaP)
	% solve the Poisson equation:   d^2 solution / dx^2 + d^2 solution/dy^2 = rhs
	% Numerically, this equation results in a matrix inversion:
	% A_ij solution_j = rhs_i
	% for more details, see Section 3.4 on page 8
	%
	% inputs:
	% rhs: right hand side to Poisson equation, size [N+2*ng   N+2*ng]
	% h: grid size
	% BC: array defining boundary conditions, size [4]
	% ng: number of guardcells
	% N: number of points in the grid
	% A: matrix defining Poisson operator
	%
	% returns:
	% solution: the solution to the Poisson equation, the new pressure, size [N+2*ng   N+2*ng]
	solution = rhs;
	
	% make the vector f that defines our right hand side	
	%
	% this is equivalent to
	%bigM = N^2;
	%for i=1:N; for j=1:N
	%	m(i,j) = i + (j-1)*N;
	%	Mlookup( m(i,j), :) = [i,j];	
	%end; end
	%for Midx = 1:bigM
	%	idx = Mlookup(Midx,:);
	%	i = idx(1); j = idx(2);
	%	f( Midx,1 ) = rhs(i+ng,j+ng); 
        %end
        f = reshape(rhs(ng+1:ng+N,ng+1:ng+N),N^2,1);   
                                	           
	% solve the Poisson equation for pressure Ap = f
	if length(LaP.perp)>0
		% if we can do the Cholesky decomposition, solve that way
        	solnVector(LaP.perp) = -LaP.Rp\(LaP.Rp'\f(LaP.perp));
	else
		% if not, use matlab's basic linear solve
		solnVector = -LaP.A\f;
	end
                                	           
	% convert back to x,y,z coordinates
	%
	% this is equivalent to
	%for Midx = 1:bigM
	%	idx = Mlookup(Midx,:);
	%	i = idx(1); j=idx(2);
	%	solution(i+ng,j+ng) = solnVector(Midx);
	%end
        solution(ng+1:ng+N,ng+1:ng+N) = reshape(solnVector,N,N);

end
%--------------------------------------------------------------------------------------------------
% functions called by poissonSolve:
% make_matrix (in own file, make_matrix.m)
% \ (Matlab function)
%
%--------------------------------------------------------------------------------------------------
