function [A,f,pGuess] = make_matrix(bigM,N1,N2,ng,h,Mlookup,m,rhs,p,BC)
	% form the matrix inversion problem that solves for the pressure
	% Ap = f
	A = sparse(bigM,bigM);
	for Midx = 1:bigM
		A(Midx,Midx) = -4;
		idx = Mlookup(Midx,:);
		i = idx(1); j = idx(2);
		f( m(i,j),1 ) = (h^2) * rhs(i+ng,j+ng); 
		pGuess( m(i,j),1 ) = p(i+ng,j+ng);
		% interior point---------------------------------------------------------------------------------
		if i>1 && i < N1 &&  j>1 && j < N2 
			A(Midx, m(i-1,j) ) = 1; 
			A(Midx, m(i+1,j) ) = 1; 
			A(Midx, m(i,j-1) ) = 1; 
			A(Midx, m(i,j+1) ) = 1; 
		% lower x boundary-------------------------------------------------------------------------------
		elseif  i ==1 &&  j>1 && j < N2 
			if BC(1)== sqrt(-1) % periodic BC
				A(Midx, m(N1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(2,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 	
			endif		
		% upper x boundary-------------------------------------------------------------------------------
		elseif i==N1 &&  j>1 && j < N2
			if BC(2)==sqrt(-1) % periodic BC
                                        			A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(N1-1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			endif
		% lower y boundary-------------------------------------------------------------------------------
		elseif i>1 && i < N1 && j==1 
			if BC(3)==sqrt(-1) % periodic BC
                                        			A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,N2) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else  % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,2) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			endif
		% upper y boundary-------------------------------------------------------------------------------
		elseif i>1 && i < N1 &&  j==N2 
			if BC(4)==sqrt(-1) % periodic BC
                                        			A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,1) ) = 1; 
			else  % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,N2-1) ) = 1; 
			endif
		% lower left corner------------------------------------------------------------------------------
		elseif i==1 && j==1 
			if BC(1)==sqrt(-1) % periodic  BC
				%A(Midx, m(N1,j) ) = 1; 
				%A(Midx, m(i+1,j) ) = 1; 
				%A(Midx, m(i,N2) ) = 1; 
				%A(Midx, m(i,j+1) ) = 1; 
				% pin down corner to allow matrix to solve
				A(Midx,Midx) = 1;
				f(Midx,1) = 1;
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(2,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,2) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			endif
		% upper left corner------------------------------------------------------------------------------
		elseif i==1 && j==N2 
			if BC(1)==sqrt(-1) % periodic BC
				A(Midx, m(N1,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(2,j) ) = 1; 
				A(Midx, m(i+1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,N2-1) ) = 1; 
			endif
		% lower right corner-----------------------------------------------------------------------------
		elseif i==N1 && j==1 
			if BC(1)==sqrt(-1) % periodic BC
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(1,j) ) = 1; 
				A(Midx, m(i,N2) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(N1-1,j) ) = 1; 
				A(Midx, m(i,2) ) = 1; 
				A(Midx, m(i,j+1) ) = 1; 
			endif
		% upper right corner-----------------------------------------------------------------------------
		elseif i==N1 && j==N2 
			if BC(1)==sqrt(-1) % periodic BC
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,1) ) = 1; 
			else % BC is wall moving with speed = BC parallel to itself
				A(Midx, m(i-1,j) ) = 1; 
				A(Midx, m(N1-1,j) ) = 1; 
				A(Midx, m(i,j-1) ) = 1; 
				A(Midx, m(i,N2-1) ) = 1; 
			endif
		endif
	endfor
	
endfunction