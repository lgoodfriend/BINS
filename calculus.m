% calculus
% this file contains routines that estimate derivatives:
%	adv_duudx
%	adv_duvdy
%	adv_dvudx
%	adv_dvvdy
%	firstDeriv
%	secondDeriv
%	divergence
%	correct_dpdx
%	correct_dpdy

1; % this is needed in octave for the interpreter to find all the functions
%------------------------------------------------------------------------------------------------------------------------------------------
% all of the advection term functions
% these are very similar, but all differ slightly due to the staggered grid
% for an overview of the advection term derivatives, see Section 3.3.2 on page 6
function duudx = adv_duudx(u,h)
	% calculate the quantity d(uu)/dx for use in the x-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% u: x-direction velocity (N+2*ng+1 x N+2*ng array)
	% h: spatial step in x direction
	%
	% returns:
	% duudx: derivative of quantity uu with respect to x (N+2*ng+1 x N+2*ng array)
	
	% initialize to get correct shape
	duudx = zeros(size(u));
	% average u to cell centers
	u_right = 0.5*(u(3:end,:)+u(2:end-1,:));
	u_left = 0.5*(u(2:end-1,:)+u(1:end-2,:));
	% combine
	uu_right = u_right.*u_right;
	uu_left = u_left.*u_left;
	% take derivative
	duudx(2:end-1,:) = (uu_right-uu_left)./h;
endfunction
function duvdy = adv_duvdy(u,v,h)
	% calculate the quantity d(uv)/dy for use in the x-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% u: x-direction velocity (N+2*ng+1 x N+2*ng array)
	% v: y-direction velocity (N+2*ng x N+2*ng+1 array)
	% h: spatial step in y direction
	%
	% returns:
	% duvdy: derivative of quantity uv with respect to y (N+2*ng+1 x N+2*ng array)
	
	% initialize to get correct shape
	duvdy = zeros(size(u));
	% average u to corners
	u_up = 0.5*(u(2:end-1,3:end)+u(2:end-1,2:end-1));
	u_down = 0.5*(u(2:end-1,2:end-1)+u(2:end-1,1:end-2));
	% average v to corners
	v_up = 0.5*(v(1:end-1,3:end-1)+v(2:end,3:end-1));
	v_down = 0.5*(v(1:end-1,2:end-2)+v(2:end,2:end-2));
	% combine
	uv_up = u_up.*v_up;
	uv_down = u_down.*v_down;
	% take derivative
	duvdy(2:end-1,2:end-1) = (uv_up-uv_down)/h;
endfunction
function dvudx = adv_dvudx(u,v,h)
	% calculate the quantity d(vu)/dx for use in the y-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% u: x-direction velocity (N+2*ng+1 x N+2*ng array)
	% v: y-direction velocity (N+2*ng x N+2*ng+1 array)
	% h: spatial step in y direction
	%
	% returns:
	% dvudx: derivative of quantity vu with respect to y  (N+2*ng x N+2*ng+1 array)
	
	% initialize to get correct shape
	dvudx = zeros(size(v));
	% average v to corners
	v_right = 0.5*(v(3:end,2:end-1)+v(2:end-1,2:end-1));
	v_left = 0.5*(v(2:end-1,2:end-1)+v(1:end-2,2:end-1));
	% average u to corners
	u_right = 0.5*(u(3:end-1,2:end)+u(3:end-1,1:end-1));
	u_left = 0.5*(u(2:end-2,2:end)+u(2:end-2,1:end-1));
	% combine
	vu_right = v_right.*u_right;
	vu_left = v_left.*u_left;
	% take derivative
	dvudx(2:end-1,2:end-1) = (vu_right-vu_left)/h;
endfunction
function dvvdy = adv_dvvdy(v,h)
	% calculate the quantity d(vv)/dy for use in the y-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% v: y-direction velocity  (N+2*ng x N+2*ng+1 array)
	% h: spatial step in y direction
	%
	% returns:
	% dvvdy: derivative of quantity vv with respect to y  (N+2*ng x N+2*ng+1 array)
	
	% initialize to get correct shape
	dvvdy = zeros(size(v));
	% average v to cell centers
	v_up = 0.5*(v(:,3:end)+v(:,2:end-1));
	v_down = 0.5*(v(:,2:end-1)+v(:,1:end-2));
	% combine
	vv_up = v_up.*v_up;
	vv_down = v_down.*v_down;
	% take derivative
	dvvdy(:,2:end-1) = (vv_up-vv_down)./h;
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
function derivOut = secondDeriv(u,dir,h)
	% calculate second derivative
	% This is much less complicated than the firstDeriv function because all output derivatives are 
	% located in the same part of the cell as the input variable.
	% for details, see Section 3.3.2 (page 6)
	%
	% inputs:
	% u: variable to take derivative of (variable dimension)
	% dir: direction of derivative, 1 for x direction, 2 for y direction
	% h: spatial step 
	%
	% returns:
	% derivOut: second derivative of input variable u (same size as u)

	derivOut = u;
	switch dir
		case{1} % x direction
			derivOut(2:end-1,:) = (u(1:end-2,:) + u(3:end,:) - 2*u(2:end-1,:))/(h^2);
		case{2} % y direction
			derivOut(:,2:end-1) = (u(:,1:end-2) + u(:,3:end) - 2*u(:,2:end-1))/(h^2);
	endswitch % dir
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function divOut = divergence(u,v,h)
	% calculates the divergence
	% for an overview of first derivatives, see Section 3.3.2 on page 6
	% for details, see Section 4.1 on page 10
	%
	% inputs:
	% u: x-direction velocity (N+2*ng+1 x N+2*ng array)
	% v: y-direction velocity (N+2*ng x N+2*ng+1 array)
	% h: spatial step
	%
	% returns:
	% divOut: divergence of u and v (N+2*ng x N+2*ng array)
	
	dudx = (u(2:end,:)-u(1:end-1,:))/h;
	dvdy = (v(:,2:end)-v(:,1:end-1))/h;
	divOut = dudx + dvdy;
endfunction
%------------------------------------------------------------------------------------------------------------------------------------------
% both of the pressure derivative functions
% like the advection terms, these are very similar but differ due to the staggered grid
% for an overview of the first derivatives, see Section 3.3.2 on page 6
function dpdx = correct_dpdx(p,h)
	% calculates dpdx on cell faces from p on cell centers
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% p: pressure (N+2*ng x N+2*ng array)
	% h: spatial step
	%
	% returns:
	% dpdx: derivative of pressure in the x direction (N+2*ng+1 x N+2*ng array)
	
	% initialize to correct size
	size_p = size(p);
	dpdx = zeros( size_p(1)+1, size_p(2));
	% take derivative
	dpdx(2:end-1,:) = (p(2:end,:)-p(1:end-1,:))/h;
endfunction
function dpdy = correct_dpdy(p,h)
	% calculates dpdx on cell faces from p on cell centers
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% p: pressure (N+2*ng x N+2*ng array)
	% h: spatial step
	%
	% returns:
	% dpdy: derivative of pressure in the y direction (N x N+1 array)
	
	% initialize to correct size
	size_p = size(p);
	dpdy = zeros( size_p(1), size_p(2)+1);
	% take derivative
	dpdy(:,2:end-1) = (p(:,2:end)-p(:,1:end-1))/h;
endfunction




