function [u,v,p] = IC(L,N,ng,IC_choice)
	% base size of variables is N+2*ng
	% each dimension looks like
	% [guardcells     domain        guardcells]
	%
	% staggered grid means u and v are offset:
	%  ______________
	% |                        |
	% |                        |
	% |                        |
	% |->u      p          |
	% |                        |
	% |           ^ v       |
	% |_______|______|
	%
	u = zeros(N(1)+2*ng+1, N(2)+2*ng);
	v = zeros(N(1)+2*ng, N(2)+2*ng+1);
	p = zeros(N(1)+2*ng, N(2)+2*ng);
	
	switch IC_choice
		case{1} % all zeros
		% no change, initialization is all zeros
		
		case{2} % velocity is all ones
		u = u+1; v = v+1;
		
		case{3} % Taylor-Green vortex
		h = L./(N-1);
		for i=ng+1:ng+N(1)+1; for j=ng+1:ng+N(2)+1
			% determine (x,y) coordinate of (i,j) point: centered
			x = (i-ng-1)*h(1); y = (j-ng-1)*h(2);
			p(i,j) = (1/4)*( cos(2*x) + cos(2*y) );
			% determine (x,y) coordinate of (i,j) point: left face
			x = (i-ng-1-0.5)*h(1); y = (j-ng-1)*h(2);
			u(i,j) = sin(x)*cos(y);
			% determine (x,y) coordinate of (i,j) point: bottom face
			x = (i-ng-1)*h(1); y = (j-ng-1-0.5)*h(2);
			v(i,j) = -cos(x)*sin(y);
		endfor; endfor;
		
	endswitch
	
endfunction