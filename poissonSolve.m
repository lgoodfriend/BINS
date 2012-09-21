function solution = poissonSolve(rhs,p,h,BC,ng)
	% solve the Poisson equation d^2 solution / dx_i dx_i = rhs
	%
	% This equation results in a matrix inversion:
	% A_ij soln_j = rhs_i

	% find the dimensions of the matrix
	N = size(rhs);
	M = prod(N);
	solution = p;

	N1 = N(1) - 2*ng; N2 = N(2) - 2*ng;
	bigM = N1*N2;
	A =sparse(bigM,bigM); f = sparse(bigM,1); pGuess = zeros(bigM,1);
	for i=1:N1; for j=1:N2
		m(i,j) = i + (j-1)*N1;
		Mlookup( m(i,j), :) = [i,j];                    	          	
	endfor; endfor
          	                          	
	[A,f,pGuess] = make_matrix(bigM,N1,N2,ng,h,Mlookup,m,rhs,p,BC);
                                	           
	solnVector = A\f;
                                	           
	% convert back to x,y,z coordinates
	for Midx = 1:bigM
		idx = Mlookup(Midx,:);
		i = idx(1); j=idx(2);
		solution(i+ng,j+ng) = solnVector(Midx);
	endfor
endfunction







