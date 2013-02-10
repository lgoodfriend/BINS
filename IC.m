function [u,v,p] = IC(L,N,ng,IC_choice)
	% initialize u,v, and p
	% initial conditions are discussed in section 2.3 on page 4
	%
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
	% u: x-direction velocity, size [N+2*ng+1   N+2*ng]
	% v: y-direction velocity, size [N+2*ng   N+2*ng+1]
	% p: pressure, size [N+2*ng   N+2*ng]
	
	u = zeros(N+2*ng+1, N+2*ng);
	v = zeros(N+2*ng, N+2*ng+1);
	p = zeros(N+2*ng, N+2*ng);
	
	h = L/N;
	li = ng+1; % lower index of inner (non-ghost cell) data
	ui = ng+N; % upper index of inner (non-ghost cell) data
	
	switch IC_choice
		case{1} % all zeros
		% no change, initialization is all zeros
		
		case{2} % velocity is all ones
		u = u+1; v = v+1;
		
		case{3} % Taylor-Green vortex
		% u
		xx = linspace(0,L,N+1);
		yy = linspace(h/2,L-h/2,N);
		[y,x] = meshgrid(yy,xx);
		u(li:ui+1,li:ui) = sin(x).*cos(y);
		% v
		xx = linspace(h/2,L-h/2,N);
		yy = linspace(0,L,N+1);
		[y,x] = meshgrid(yy,xx);
		v(li:ui,li:ui+1) = -cos(x).*sin(y);
		% p
		xx = linspace(h/2,L-h/2,N);
		yy = linspace(h/2,L-h/2,N);
		[y,x] = meshgrid(yy,xx);
		p(li:ui,li:ui) = 0.25*(cos(2*x) + cos(2*y));
		
		case{4} % shear flow
		% u
		xx = linspace(0,L,N+1);
		yy = linspace(h/2,L-h/2,N);
		[y,x] = meshgrid(yy,xx);
		u(li:ui+1,li:li+N/2-1) = tanh(30*(y(:,1:N/2)-1/4));
		u(li:ui+1,li+N/2:ui) = tanh(30*(3/4-y(:,N/2+1:N)));
		% v
		xx = linspace(h/2,L-h/2,N);
		yy = linspace(0,L,N+1);
		[y,x] = meshgrid(yy,xx);
		v(li:ui,li:ui+1) = (1/20)*sin(2*pi*x);
		% p
		p=zeros(size(p));
		
	endswitch
	
endfunction