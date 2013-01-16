function [u,v,p] = IC(L,N,ng,IC_choice)
	% base size of variables is N+2*ng
	% each dimension looks like
	% [guardcells     domain        guardcells]
	%
	% staggered grid means u and v are offset:
	%  ______________
	% |                             |
	% |                             |
	% |                             |
	% |->u      p              |
	% |                             |
	% |              ^ v         |
	% |_______|______|
	%
	% inputs:
	% L: length of domain
	% N: number of grid cells
	% ng: number of guardcells
	% IC_choice: integer to choose among initial conditions defined here
	%
	% returns:
	% u: x-direction velocity (N+2*ng+1 x N+2*ng array)
	% v: y-direction velocity (N+2*ng x N+2*ng+1 array)
	% p: pressure (N+2*ng x N+2*ng array)
	
	u = zeros(N+2*ng+1, N+2*ng);
	v = zeros(N+2*ng, N+2*ng+1);
	p = zeros(N+2*ng, N+2*ng);
	
	switch IC_choice
		case{1} % all zeros
		% no change, initialization is all zeros
		
		case{2} % velocity is all ones
		u = u+1; v = v+1;
		
		case{3} % Taylor-Green vortex
		h = L./(N-1);
		for i=ng+1:ng+N+1; for j=ng+1:ng+N+1
			% determine (x,y) coordinate of (i,j) point: centered
			x = (i-ng-1 + 0.5)*h; y = (j-ng-1 + 0.5)*h;
			p(i,j) = (1/4)*( cos(2*x) + cos(2*y) );
			% determine (x,y) coordinate of (i,j) point: left face
			x = (i-ng-1)*h; y = (j-ng-1 + 0.5)*h;
			u(i,j) = sin(x)*cos(y);
			% determine (x,y) coordinate of (i,j) point: bottom face
			x = (i-ng-1 + 0.5)*h; y = (j-ng-1)*h;
			v(i,j) = -cos(x)*sin(y);
		endfor; endfor;
		
		case{4} % shear flow
		h = L./(N-1);
		for i=ng+1:ng+N+1; for j=ng+1:ng+N+1
			% determine (x,y) coordinate of (i,j) point: centered
			x = (i-ng-1 + 0.5)*h; y = (j-ng-1 + 0.5)*h;
			p(i,j) = 0;
			% determine (x,y) coordinate of (i,j) point: left face
			x = (i-ng-1)*h; y = (j-ng-1 + 0.5)*h;
			if y<=L/2 
				u(i,j) = tanh(30*(y-1/4));
			else
				u(i,j) = tanh(30*(3/4-y));
			endif
			% determine (x,y) coordinate of (i,j) point: bottom face
			x = (i-ng-1 + 0.5)*h; y = (j-ng-1)*h;
			v(i,j) = (1/20)*sin(2*pi*x);
		endfor; endfor;
		
	endswitch
	
endfunction