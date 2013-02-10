% advance
% this file contains functions that advance the solution in time:
%	advection
%	diffusion
%	predict
%	pressure
%	correct

1; % this is needed in octave for the interpreter to find all the functions
%------------------------------------------------------------------------------------------------------------------------------------------
function [u_advection, v_advection] = advection(u,v,h)
% calculates the advection term of the Navier-Stokes equation
% the advection term is discussed in Section 3.2, starting on page 4
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

endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [u_diffusion,v_diffusion] = diffusion(u,v,h,nu)
% calculates the diffusion term of the Navier-Stokes equation
% the diffusion term is discussed in Section 3.2, starting on page 4
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
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [uStar,vStar] = predict(u,v,u_advection,v_advection,u_diffusion,v_diffusion,dt)
% estimates the velocity at the next time step
% the prediction step is discussed in Section 3.2, starting on page 4
% the time stepping method used here (forward Euler) is discussed in Section 3.3.1 on page 6
% this function solves equations 1 and 2 for u* (page 5)
%
% inputs:
% u = x-direction velocity, size [N+2*ng+1   N+2*ng]
% v = y-direction velocity, size [N+2*ng   N+2*ng+1]
% u_advection = x-direction advection term, size [N+2*ng+1   N+2*ng]
% v_advection = y-direction advection term, size [N+2*ng   N+2*ng+1]
% u_diffusion = x-direction diffusion term, size [N+2*ng+1   N+2*ng]
% v_diffusion = y-direction diffusion term, size [N+2*ng   N+2*ng+1]
% dt = timestep
%
% returns:
% uStar = x-direction predicted velocity at next time step, size [N+2*ng+1   N+2*ng]
% vStar = y-direction predicted velocity at next time step, size [N+2*ng   N+2*ng+1]
	uStar = u + (-u_advection + u_diffusion)*dt;
	vStar = v + (-v_advection + v_diffusion)*dt;
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function newP = pressure(uStar,vStar,p,h,dt,ng,BC,N)
% solves for the pressure using the requirement that the velocity
% field must be divergence free:
% laplacian(p) =  divergence(uStar) / dt
% the Poisson pressure solve is discussed in Section 3.2, starting on page 4 
% this function solves equation 3 for p^n+1 (page 5)
%
% inputs:
% uStar = predicted x-direction velocity at next time step, size [N+2*ng+1   N+2*ng]
% vStar = predicted y-direction velocity at next time step, size [N+2*ng   N+2*ng+1]
% p = pressure, size [N+2*ng  N+2*ng]
% h = dx and dy
% dt = timestep
% ng = number of guardcells
% BC = boundary conditions, size [4]
% N = number of cells in each direction
%
% returns:
% newP = pressure at next time step
	newP = p;
	% calculate right hand side of Poisson equation
	div = divergence(uStar,vStar,h);
	rhs = div/dt;
	% solve Poisson equation for pressure
	newP = poissonSolve(rhs,h,BC,ng,N);
	
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [uNext,vNext] = correct(uStar,vStar,newP,h,dt)
% correct the velocity with the new pressure field 
% this forces the velocity field to be divergence free
% the correction step is described in Section 3.2, starting on page 4
% this function solves equations 4 and 5 (page 5)
%
% inputs:
% uStar = predicted x-direction velocity at next time step, size [N+2*ng+1   N+2*ng]
% vStar = predicted y-direction velocity at next time step, size [N+2*ng   N+2*ng+1]
% newP = pressure at next time step
% h = dx and dy
% dt = timestep
%
% returns:
% uNext = x-direction velocity at next time step, size [N+2*ng+1   N+2*ng]
% vNext = y-direction velocity at next time step, size [N+2*ng   N+2*ng+1]
	% calculate new pressure derivatives
	dpdx = correct_dpdx(newP,h);
	dpdy = correct_dpdy(newP,h);
	% correct velocity to be divergence-free with new pressure derivatives
	uNext = uStar - dpdx*dt;
	vNext = vStar - dpdy*dt;

endfunction
