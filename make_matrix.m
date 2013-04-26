function A = make_matrix(N,ng,h,BC)
	% form the matrix inversion problem that solves for the pressure
	% Ap = f
	% for details, see Section 3.4 on page 8
	%
	% inputs:
	% bigM: the size of the matrix
	% N: number of grid points in x and y direction
	% ng: number of guardcells
	% h: grid spacing dx and dy
	% Mlookup: matrix to look up (i,j) coordinates from matrix location, size [bigM 2]
	% m: matrix to look up matrix location from (i,j) coordinates, size [N  N]
	% rhs: right hand side of Poisson equation, size [N+2*ng   N+2*ng]
	% BC: array of boundary conditions, size [4]
	%
	% returns:
	% A: matrix defining Poisson operator and boundary conditions, size [bigM   bigM]
	% f: vector defining right hand side, size [bigM]

	% make a lookup table for converting between matrix and vector representations
	% of the xy plane
	bigM = N^2;
	for i=1:N; for j=1:N
		m(i,j) = i + (j-1)*N;
		Mlookup( m(i,j), :) = [i,j];                    	          	
	end; end
	
	A = sparse(bigM,bigM);
	f = sparse(bigM,1);
	for Midx = 1:bigM
		A(Midx,Midx) = -4;
		idx = Mlookup(Midx,:);
		i = idx(1); j = idx(2);
		% interior point---------------------------------------------------------------------------------
		if i>1 && i < N &&  j>1 && j < N 
			A(Midx, m(i-1,j) ) = 1; 
			A(Midx, m(i+1,j) ) = 1; 
			A(Midx, m(i,j-1) ) = 1; 
			A(Midx, m(i,j+1) ) = 1; 
		% -----------------------------------------------------------------------------------------------
		% boundary conditions----------------------------------------------------------------------------
		% -----------------------------------------------------------------------------------------------
		% lower x boundary-------------------------------------------------------------------------------
		elseif  i ==1 &&  j>1 && j < N 
			if BC(1)== sqrt(-1) % periodic BC
				A(Midx, m(N,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 	
				A(Midx,Midx) = A(Midx,Midx)+1;
			end	
		% upper x boundary-------------------------------------------------------------------------------
		elseif i==N &&  j>1 && j < N
			if BC(2)==sqrt(-1) % periodic BC
                                A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
				A(Midx,Midx) = A(Midx,Midx)+1;
			end
		% lower y boundary-------------------------------------------------------------------------------
		elseif i>1 && i < N && j==1 
			if BC(3)==sqrt(-1) % periodic BC
                                A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,N) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else  % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
				A(Midx,Midx) = A(Midx,Midx)+1;
			end
		% upper y boundary-------------------------------------------------------------------------------
		elseif i>1 && i < N &&  j==N 
			if BC(4)==sqrt(-1) % periodic BC
                                A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,1) ) = 1; 
			else  % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, Midx) = A(Midx,Midx)+1;
			end
		% lower left corner------------------------------------------------------------------------------
		elseif i==1 && j==1 
                        % this is what the BCs would be, except that we must "pin down" one corner of the matrix
                        % to make it invertible when using periodic or wall BCs
			%% x BCs
			%if BC(1)==sqrt(-1) % periodic  BC in x
		        %	A(Midx, m(N,j) ) = 1; 
			%	A(Midx, m(i+1,j) ) = 1; 
			%else % x BC is wall moving with speed = BC parallel to itself
			%	A(Midx, Midx ) = A(Midx,Midx)+1; 
			%	A(Midx, m(i+1,j) ) = 1; 
			%end
			%% y BCs
			%if BC(3)==sqrt(-1) % periodic BC in y
			%	A(Midx, m(i,N) ) = 1;
			%	A(Midx, m(i,j+1) ) = 1;
			%else % y BC is wall moving with speed = BC parallel to itself
			%	A(Midx, Midx) = A(Midx,Midx)+1;
			%	A(Midx, m(i,j+1) ) = 1;
			%end
                        A(Midx,Midx) = 1;
                        f(Midx) = 1;  
		% upper left corner------------------------------------------------------------------------------
		elseif i==1 && j==N 
			% x BCs
			if BC(1)==sqrt(-1) % periodic  BC in x
				A(Midx, m(N,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
			else % x BC is wall moving with speed = BC parallel to itself
				A(Midx, Midx ) = A(Midx,Midx)+1; 
				A(Midx, m(i+1,j) ) = 1; 
			end
			% y BCs
			if BC(3)==sqrt(-1) % periodic BC in y
				A(Midx, m(i,1) ) = 1;
				A(Midx, m(i,j-1) ) = 1;
			else % y BC is wall moving with speed = BC parallel to itself
				A(Midx, Midx) = A(Midx,Midx)+1;
				A(Midx, m(i,j-1) ) = 1;
			end
		% lower right corner-----------------------------------------------------------------------------
		elseif i==N && j==1 
			% x BCs
			if BC(1)==sqrt(-1) % periodic  BC in x
				A(Midx, m(1,j) ) = 1; 
				A(Midx, m(i-1,j) ) = 1; 
			else % x BC is wall moving with speed = BC parallel to itself
				A(Midx, Midx ) = A(Midx,Midx)+1; 
				A(Midx, m(i-1,j) ) = 1; 
			end
			% y BCs
			if BC(3)==sqrt(-1) % periodic BC in y
				A(Midx, m(i,N) ) = 1;
				A(Midx, m(i,j+1) ) = 1;
			else % y BC is wall moving with speed = BC parallel to itself
				A(Midx, Midx) = A(Midx,Midx)+1;
				A(Midx, m(i,j+1) ) = 1;
			end
		% upper right corner-----------------------------------------------------------------------------
		elseif i==N && j==N 
			% x BCs
			if BC(1)==sqrt(-1) % periodic  BC in x
				A(Midx, m(1,j) ) = 1; 
				A(Midx, m(i-1,j) ) = 1; 
			else % x BC is wall moving with speed = BC parallel to itself
				A(Midx, Midx ) = A(Midx,Midx)+1; 
				A(Midx, m(i-1,j) ) = 1; 
			end
			% y BCs
			if BC(3)==sqrt(-1) % periodic BC in y
				A(Midx, m(i,1) ) = 1;
				A(Midx, m(i,j-1) ) = 1;
			else % y BC is wall moving with speed = BC parallel to itself
				A(Midx, Midx) = A(Midx,Midx)+1;
				A(Midx, m(i,j-1) ) = 1;
			end
		end
	end % loop over rows in matrix A
	
end
