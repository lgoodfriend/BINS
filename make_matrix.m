function poisSolve = make_matrix(N,h,BC)
	% form the matrix inversion problem that solves for the pressure
	% Laplacian(p) = f
	% for details, see Section 3.4 on page 8
	%
	% inputs:
	% N: number of grid points in x and y direction
	% h: grid spacing dx and dy
	% BC: array of boundary conditions, size [4]
	%
	% returns:
	% poisSolve: a structure containing either
	%	if the BCs make the Laplacian operator positive definite 
	%	(true for the wall BCs)
	%	perp: a sparse matrix defining how the indices of Rp have been rearranged
	%	Rp: a sparse matrix defining the Laplace operator with rearranged indices
	%
	%	or, if the BCs do not,
	%	(true for the periodic BCs)
	%	A: a sparse matrix defining the Laplace operator

        % form the matrix defining the Laplacian operator
	% see below for the slow but clear way of forming this matrix 
        A = kron(speye(N), LaP1D(N,h,BC(1)) ) + kron( LaP1D(N,h,BC(3)) ,speye(N));

	% reorder the terms for a fast Cholesky decomposition
        poisSolve.perp = symamd(A); 
        [poisSolve.Rp,err] = chol(A(poisSolve.perp,poisSolve.perp));

        if err>0 % then A wasn't positive definite, so can't do a Cholesky decomposition
		poisSolve.perp = [];
		poisSolve.Rp = [];
		A(1,:) = sparse(1,N^2); A(1,1) = -1/h^2; % pin corner so matrix is solvable
		poisSolve.A = A;
	end

        % The code above is just a compact way of building the Laplacian operator matrix
        % with the correct boundary conditions. It is equivalent to the following,
        % which is both more lines of code and slower (in Matlab).
        %
	% make a lookup table for converting between matrix and vector representations
	% of the xy plane
%	bigM = N^2;
%	for i=1:N; for j=1:N
%		m(i,j) = i + (j-1)*N;
%		Mlookup( m(i,j), :) = [i,j];                    	          	
%	end; end
%	
%	A = sparse(bigM,bigM);
%	f = sparse(bigM,1);
%	for Midx = 1:bigM
%		A(Midx,Midx) = -4;
%		idx = Mlookup(Midx,:);
%		i = idx(1); j = idx(2);
%		% interior point-------------------------------------------------------------
%		if i>1 && i < N &&  j>1 && j < N 
%			A(Midx, m(i-1,j) ) = 1; 
%			A(Midx, m(i+1,j) ) = 1; 
%			A(Midx, m(i,j-1) ) = 1; 
%			A(Midx, m(i,j+1) ) = 1; 
%		% ----------------------------------------------------------------------------
%		% boundary conditions---------------------------------------------------------
%		% ----------------------------------------------------------------------------
%		% lower x boundary------------------------------------------------------------
%		elseif  i ==1 &&  j>1 && j < N 
%			if BC(1)== sqrt(-1) % periodic BC
%				A(Midx, m(N,j) ) = 1; 
%				A(Midx, m(i+1,j) ) = 1; 
%				A(Midx, m(i,j-1) ) = 1; 
%				A(Midx, m(i,j+1) ) = 1; 
%			else % BC is wall moving with speed = BC parallel to itself
%				A(Midx, m(i+1,j) ) = 1; 
%				A(Midx, m(i,j-1) ) = 1; 
%				A(Midx, m(i,j+1) ) = 1; 	
%				A(Midx,Midx) = A(Midx,Midx)+1;
%			end	
%		% upper x boundary-------------------------------------------------------------
%		elseif i==N &&  j>1 && j < N
%			if BC(2)==sqrt(-1) % periodic BC
%                                A(Midx, m(i-1,j) ) = 1; 
%				A(Midx, m(1,j) ) = 1; 
%				A(Midx, m(i,j-1) ) = 1; 
%				A(Midx, m(i,j+1) ) = 1; 
%			else % BC is wall moving with speed = BC parallel to itself
%				A(Midx, m(i-1,j) ) = 1; 
%				A(Midx, m(i,j-1) ) = 1; 
%				A(Midx, m(i,j+1) ) = 1; 
%				A(Midx,Midx) = A(Midx,Midx)+1;
%			end
%		% lower y boundary-------------------------------------------------------------
%		elseif i>1 && i < N && j==1 
%			if BC(3)==sqrt(-1) % periodic BC
%                                A(Midx, m(i-1,j) ) = 1; 
%				A(Midx, m(i+1,j) ) = 1; 
%				A(Midx, m(i,N) ) = 1; 
%				A(Midx, m(i,j+1) ) = 1; 
%			else  % BC is wall moving with speed = BC parallel to itself
%				A(Midx, m(i-1,j) ) = 1; 
%				A(Midx, m(i+1,j) ) = 1; 
%				A(Midx, m(i,j+1) ) = 1; 
%				A(Midx,Midx) = A(Midx,Midx)+1;
%			end
%		% upper y boundary-------------------------------------------------------------
%		elseif i>1 && i < N &&  j==N 
%			if BC(4)==sqrt(-1) % periodic BC
%                                A(Midx, m(i-1,j) ) = 1; 
%				A(Midx, m(i+1,j) ) = 1; 
%				A(Midx, m(i,j-1) ) = 1; 
%				A(Midx, m(i,1) ) = 1; 
%			else  % BC is wall moving with speed = BC parallel to itself
%				A(Midx, m(i-1,j) ) = 1; 
%				A(Midx, m(i+1,j) ) = 1; 
%				A(Midx, m(i,j-1) ) = 1; 
%				A(Midx, Midx) = A(Midx,Midx)+1;
%			end
%		% lower left corner------------------------------------------------------------
%		elseif i==1 && j==1 
%                        % this is what the BCs would be, except that we must 
%                        % "pin down" one corner of the matrix
%                        % to make it invertible when using periodic or wall BCs
%			%% x BCs
%			%if BC(1)==sqrt(-1) % periodic  BC in x
%		        %	A(Midx, m(N,j) ) = 1; 
%			%	A(Midx, m(i+1,j) ) = 1; 
%			%else % x BC is wall moving with speed = BC parallel to itself
%			%	A(Midx, Midx ) = A(Midx,Midx)+1; 
%			%	A(Midx, m(i+1,j) ) = 1; 
%			%end
%			%% y BCs
%			%if BC(3)==sqrt(-1) % periodic BC in y
%			%	A(Midx, m(i,N) ) = 1;
%			%	A(Midx, m(i,j+1) ) = 1;
%			%else % y BC is wall moving with speed = BC parallel to itself
%			%	A(Midx, Midx) = A(Midx,Midx)+1;
%			%	A(Midx, m(i,j+1) ) = 1;
%			%end
%                        A(Midx,Midx) = 1;
%                        f(Midx) = 1;  
%		% upper left corner----------------------------------------------------------
%		elseif i==1 && j==N 
%			% x BCs
%			if BC(1)==sqrt(-1) % periodic  BC in x
%				A(Midx, m(N,j) ) = 1; 
%				A(Midx, m(i+1,j) ) = 1; 
%			else % x BC is wall moving with speed = BC parallel to itself
%				A(Midx, Midx ) = A(Midx,Midx)+1; 
%				A(Midx, m(i+1,j) ) = 1; 
%			end
%			% y BCs
%			if BC(3)==sqrt(-1) % periodic BC in y
%				A(Midx, m(i,1) ) = 1;
%				A(Midx, m(i,j-1) ) = 1;
%			else % y BC is wall moving with speed = BC parallel to itself
%				A(Midx, Midx) = A(Midx,Midx)+1;
%				A(Midx, m(i,j-1) ) = 1;
%			end
%		% lower right corner---------------------------------------------------------
%		elseif i==N && j==1 
%			% x BCs
%			if BC(1)==sqrt(-1) % periodic  BC in x
%				A(Midx, m(1,j) ) = 1; 
%				A(Midx, m(i-1,j) ) = 1; 
%			else % x BC is wall moving with speed = BC parallel to itself
%				A(Midx, Midx ) = A(Midx,Midx)+1; 
%				A(Midx, m(i-1,j) ) = 1; 
%			end
%			% y BCs
%			if BC(3)==sqrt(-1) % periodic BC in y
%				A(Midx, m(i,N) ) = 1;
%				A(Midx, m(i,j+1) ) = 1;
%			else % y BC is wall moving with speed = BC parallel to itself
%				A(Midx, Midx) = A(Midx,Midx)+1;
%				A(Midx, m(i,j+1) ) = 1;
%			end
%		% upper right corner---------------------------------------------------------
%		elseif i==N && j==N 
%			% x BCs
%			if BC(1)==sqrt(-1) % periodic  BC in x
%				A(Midx, m(1,j) ) = 1; 
%				A(Midx, m(i-1,j) ) = 1; 
%			else % x BC is wall moving with speed = BC parallel to itself
%				A(Midx, Midx ) = A(Midx,Midx)+1; 
%				A(Midx, m(i-1,j) ) = 1; 
%			end
%			% y BCs
%			if BC(3)==sqrt(-1) % periodic BC in y
%				A(Midx, m(i,1) ) = 1;
%				A(Midx, m(i,j-1) ) = 1;
%			else % y BC is wall moving with speed = BC parallel to itself
%				A(Midx, Midx) = A(Midx,Midx)+1;
%				A(Midx, m(i,j-1) ) = 1;
%			end
%		end
%	end % loop over rows in matrix A
%       A = A/h^2;
end
%--------------------------------------------------------------------------------------------
% functions called by make_matrix:
% LaP1D
%--------------------------------------------------------------------------------------------
function A = LaP1D(N,h,BC)
	% form a matrix defining the 1D Laplace operator, with BCs
	%
	% inputs:
	% N: number of grid points in x and y direction
	% h: grid spacing dx and dy
	% BC: boundary condition for this dimension
	%
	% returns:
	% A: matrix defining the 1D Laplace operator, with BCs, for this dimension,size [N N]

	if BC==sqrt(-1); % periodic BCs
		BCval = 2; % Dirichlet BC for pressure, p = (periodic p)
 		A = spdiags([-1 BCval 0;ones(N-2,1)*[-1 2 -1];0 BCval -1],-1:1,N,N)'/h^2;
		A(1,N) = -1/h^2; A(N,1) = -1/h^2; % find the periodic point
	else; % moving wall BCs 
		BCval = 1; % Neumann BC for pressure, dpdn = 0
 		A = spdiags([-1 BCval 0;ones(N-2,1)*[-1 2 -1];0 BCval -1],-1:1,N,N)'/h^2;
	end;
end







