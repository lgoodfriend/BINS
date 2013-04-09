function [u_diffusion,v_diffusion] = diffusion(u,v,h,nu)
% calculates the diffusion term of the Navier-Stokes equation
% the diffusion term is discussed in Section 3.2, starting on page 5
% this function calculates the second terms of the right hand size of equations 1 and 2 (page 5)
%
% inputs:
% u = x-direction velocity, size [N+2*ng+1   N+2*ng]
% v = y-direction velocity, size [N+2*ng   N+2*ng+1]
% h = dx and dy
% nu = kinematic viscosity
%
% returns:
% u_diffusion = x-direction diffusion term of the Navier-Stokes equation, size [N+2*ng+1   N+2*ng]
% v_diffusion = y-direction diffusion term of the Navier-Stokes equation, size [N+2*ng   N+2*ng+1]

	% calculate second derivatives of u
	du2dx2 = secondDeriv(u,1,h);
	du2dy2 = secondDeriv(u,2,h);
	% form u diffusion term
	u_diffusion = nu*(du2dx2 + du2dy2);
	
	% calculate second derivatives of v
	dv2dx2= secondDeriv(v,1,h);
	dv2dy2= secondDeriv(v,2,h);
	% form v diffusion term
	v_diffusion = nu*(dv2dx2 + dv2dy2);
end
%--------------------------------------------------------------------------------------------------
% functions called by diffusion:
% secondDeriv
%
%--------------------------------------------------------------------------------------------------
function derivOut = secondDeriv(u,dir,h)
	% calculate second derivative
	% This is much less complicated than the adv_dXXdX functions because all output derivatives are 
	% located in the same part of the cell as the input variable.
	% for details, see Section 3.3.2 (page 7)
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
	end % dir
end

