% advance
% this file contains functions that advance the solution in time:
%	advection
%	diffusion
%	predict
%	pressure
%	correct

1; % this is needed in octave for the interpreter to find all the functions
%------------------------------------------------------------------------------------------------------------------------------------------
function [u_advection, v_advection] = advection(u,v,h,ng)
% calculates the advection term of the Navier-Stokes equation
% the advection term is discussed in Section 3.2 on page 4
% this function calculates the first term of the right hand side of equation 1 (page 4)
%
% inputs:
% u = x-direction velocity (N+2*ng+1 x N+2*ng array)
% v = y-direction velocity (N+2*ng x N+2*ng+1 array)
% h = dx and dy
% ng = number of guardcells
%
% returns:
% u_advection = x-direction advection term of the Navier-Stokes equation (N+2*ng+1 x N+2*ng array)
% v_advection = y-direction advection term of the Navier-Stokes equation (N+2*ng x N+2*ng+1 array)
	
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
function [u_diffusion,v_diffusion] = diffusion(u,v,h,ng,BC,nu)
% calculates the diffusion term of the Navier-Stokes equation
% the diffusion term is discussed in Section 3.2 on page 4
% this function calculates the second term of the right hand size of equation  1 (page 4)
%
% inputs:
% u = x-direction velocity (N+2*ng+1 x N+2*ng array)
% v = y-direction velocity (N+2*ng x N+2*ng+1 array)
% h = dx and dy
% ng = number of guardcells
% BC = boundary conditions (4x1 array)
% nu = kinematic viscosity
%
% returns:
% u_diffusion = x-direction diffusion term of the Navier-Stokes equation (N+2*ng+1 x N+2*ng array)
% v_diffusion = y-direction diffusion term of the Navier-Stokes equation (N+2*ng x N+2*ng+1 array)
	du2dx2 = secondDeriv(u,1,h);
	du2dy2 = secondDeriv(u,2,h);
	u_diffusion = nu*(du2dx2 + du2dy2);
	
	dv2dx2= secondDeriv(v,1,h);
	dv2dy2= secondDeriv(v,2,h);
	v_diffusion = nu*(dv2dx2 + dv2dy2);
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [uStar,vStar] = predict(u,v,u_advection,v_advection,u_diffusion,v_diffusion,dt)
% estimates the velocity at the next time step
% the prediction step is discussed in Section 3.2 on page 4
% the time stepping method used here (forward Euler) is discussed in Section 3.3.1 on page 5
% this function solves equation 1 for u* (page 4)
%
% inputs:
% u = x-direction velocity (N+2*ng+1 x N+2*ng array)
% v = y-direction velocity (N+2*ng x N+2*ng+1 array)
% u_advection = x-direction advection term (N+2*ng+1 x N+2*ng array)
% v_advection = y-direction advection term (N+2*ng x N+2*ng+1 array)
% u_diffusion = x-direction diffusion term (N+2*ng+1 x N+2*ng array)
% v_diffusion = y-direction diffusion term (N+2*ng x N+2*ng+1 array)
% dt = timestep
%
% returns:
% uStar = x-direction predicted velocity at next time step (N+2*ng+1 x N+2*ng array)
% vStar = y-direction predicted velocity at next time step (N+2*ng x N+2*ng+1 array)
	uStar = u + (-u_advection + u_diffusion)*dt;
	vStar = v + (-v_advection + v_diffusion)*dt;
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function newP = pressure(uStar,vStar,p,h,dt,ng,BC,N)
% solves for the pressure using the requirement that the velocity
% field must be divergence free:
% laplacian(p) =  divergence(uStar) / dt
% the Poisson pressure solve is discussed in Section 3.2 on page 4 
% this function solves equation 2 for p^n+1 (page 4)
%
% inputs:
% uStar = predicted x-direction velocity at next time step (N+2*ng+1 x N+2*ng array)
% vStar = predicted y-direction velocity at next time step (N+2*ng x N+2*ng+1 array)
% p = pressure (N+2*ng x N+2*ng array)
% h = dx and dy
% dt = timestep
% ng = number of guardcells
% BC = boundary conditions (4x1 array)
% N = number of cells in each direction (2x1 array)
%
% returns:
% newP = pressure at next time step
	newP = p;
	div = divergence(uStar,vStar,h);
	rhs = div/dt;
	newP = poissonSolve(rhs,h,BC,ng,N);
	
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [uNext,vNext] = correct(uStar,vStar,newP,h,dt,ng)
% correct the velocity with the new pressure field 
% this forces the velocity field to be divergence free
% the correction step is described in Section 3.2 on page 4
% this function solves equation 3 (page 5)
%
% inputs:
% uStar = predicted x-direction velocity at next time step (N+2*ng+1 x N+2*ng array)
% vStar = predicted y-direction velocity at next time step (N+2*ng x N+2*ng+1 array)
% newP = pressure at next time step
% h = dx and dy
% dt = timestep
% ng = number of guardcells
%
% returns:
% uNext = x-direction velocity at next time step (N+2*ng+1 x N+2*ng array)
% vNext = y-direction velocity at next time step (N+2*ng x N+2*ng+1 array)
	dpdx = correct_dpdx(newP,h);
	dpdy = correct_dpdy(newP,h);
	
	uNext = uStar - dpdx*dt;
	vNext = vStar - dpdy*dt;
endfunction
