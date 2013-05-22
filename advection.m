function [u_advection, v_advection] = advection(u,v,h)
% calculates the advection term of the Navier-Stokes equation
% the advection term is discussed in Section 3.2, starting on page 5
% this function calculates the first two terms of the right hand side of equations 1 and 2 (page 5)
%
% inputs:
% u = x-direction velocity, size [N+2*ng+1   N+2*ng]
% v = y-direction velocity, size [N+2*ng   N+2*ng+1]
% h = dx and dy
%
% returns:
% u_advection = x-direction advection term of the Navier-Stokes equation, size [N+2*ng+1   N+2*ng]
% v_advection = y-direction advection term of the Navier-Stokes equation, size [N+2*ng   N+2*ng+1]
	
	% calculate derivatives
	duudx = adv_duudx(u,h);
	duvdy = adv_duvdy(u,v,h);
	dvudx = adv_dvudx(u,v,h);
	dvvdy = adv_dvvdy(v,h);
	
	% form advection term
	u_advection = duudx + duvdy;
	v_advection = dvudx + dvvdy;

end
%--------------------------------------------------------------------------------------------------
% functions called by advection:
% adv_duudx
% adv_duvdy
% adv_dvudx
% adv_dvvdy
%
% these are very similar, but differ slightly due to the staggered grid
% for an overview of the advection term derivatives, see Section 3.3.2 on page 7
%--------------------------------------------------------------------------------------------------
function duudx = adv_duudx(u,h)
	% calculate the quantity d(uu)/dx for use in the x-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% u: x-direction velocity, size [N+2*ng+1   N+2*ng]
	% h: spatial step in x direction
	%
	% returns:
	% duudx: derivative of quantity uu with respect to x, size [N+2*ng+1   N+2*ng]
	
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
end

function duvdy = adv_duvdy(u,v,h)
	% calculate the quantity d(uv)/dy for use in the x-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% u: x-direction velocity, size [N+2*ng+1   N+2*ng]
	% v: y-direction velocity, size [N+2*ng   N+2*ng+1]
	% h: spatial step in y direction
	%
	% returns:
	% duvdy: derivative of quantity uv with respect to y, size [N+2*ng+1   N+2*ng]
	
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
end

function dvudx = adv_dvudx(u,v,h)
	% calculate the quantity d(vu)/dx for use in the y-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% u: x-direction velocity, size [N+2*ng+1   N+2*ng]
	% v: y-direction velocity, size [N+2*ng   N+2*ng+1]
	% h: spatial step in y direction
	%
	% returns:
	% dvudx: derivative of quantity vu with respect to y, size  [N+2*ng   N+2*ng+1]
	
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
end

function dvvdy = adv_dvvdy(v,h)
	% calculate the quantity d(vv)/dy for use in the y-direction advection term
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% v: y-direction velocity, size  [N+2*ng   N+2*ng+1]
	% h: spatial step in y direction
	%
	% returns:
	% dvvdy: derivative of quantity vv with respect to y, size  [N+2*ng   N+2*ng+1]
	
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
end
