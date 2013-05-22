function solution = poissonSolve(rhs,h,BC,ng,N,A)
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
	% equivalent to
        % f = (h^2)*reshape(rhs(ng+1:ng+N,ng+1:ng+N),N^2,1);   
	bigM = N^2;
	for i=1:N; for j=1:N
		m(i,j) = i + (j-1)*N;
		Mlookup( m(i,j), :) = [i,j];	
	end; end
	for Midx = 1:bigM
		idx = Mlookup(Midx,:);
		i = idx(1); j = idx(2);
		f( Midx,1 ) = (h^2) * rhs(i+ng,j+ng); 
        end
                                	           
	% solve solnVector = A\f using matlab's built-ins
	solnVector = A\f;
                                	           
	% convert back to x,y,z coordinates
	%
	% equivalent to
        % solution(ng+1:ng+N,ng+1:ng+N) = reshape(solnVector,N,N);
	for Midx = 1:bigM
		idx = Mlookup(Midx,:);
		i = idx(1); j=idx(2);
		solution(i+ng,j+ng) = solnVector(Midx);
	end

end
%--------------------------------------------------------------------------------------------------
% functions called by poissonSolve:
% make_matrix (in own file, make_matrix.m)
% \ (Matlab function)
%
%--------------------------------------------------------------------------------------------------
