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
% see section XX of the documentation for more details
%
% inputs:
% u = x-direction velocity (Nx x Ny array)
% v = y-direction velocity (Nx x Ny array)
% h = dx and dy (2x1 array)
% ng = number of guardcells
%
% returns:
% advectionTerm = advection term of the Navier-Stokes equation (2 x Nx x Ny array)
	% average velocities to cell centers
	uC = 0.5*(u(1:end-1,:) + u(2:end,:));
	vC = 0.5*(v(:,1:end-1) + v(:,2:end));
	% take derivatives
	duudx = firstDeriv(uC.*uC,1,h,'advU');
	duvdy = firstDeriv(uC.*vC,2,h,'advU');
	dvudx = firstDeriv(vC.*uC,1,h,'advV');
	dvvdy = firstDeriv(vC.*vC,2,h,'advV');
	% form advection term
	u_advection = duudx + duvdy;
	v_advection = dvudx + dvvdy;
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [u_diffusion,v_diffusion] = diffusion(u,v,h,ng,BC,nu)
% calculates the diffusion term of the Navier-Stokes equation
% see section XX of the documentation for more details
%
% inputs:
% u = velocity (2 x Nx x Ny array)
% h = dx and dy (2x1 array)
% ng = number of guardcells
% BC = boundary conditions (4x1 array)
% nu = kinematic viscosity
%
% returns:
% diffusionTerm = diffusion term of the Navier-Stokes equation (2 x Nx x Ny array)
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
% see section XX of the documentation for more details
%
% inputs:
% u = velocity (2 x Nx x Ny array)
% advection = advection term (2 x Nx x Ny array)
% diffusion = diffusion term (2 x Nx x Ny array)
% dt = timestep
%
% returns:
% uStar = predicted velocity at next time step (2 x Nx x Ny array)
	uStar = u + (-u_advection + u_diffusion)*dt;
	vStar = v + (-v_advection + v_diffusion)*dt;
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function newP = pressure(uStar,vStar,p,h,dt,ng,BC,N)
% solves for the pressure using the requirement that the velocity
% field must be divergence free:
% laplacian(p) =  divergence(uStar)
% see section XX of the documentation for more details
%
% inputs:
% uStar = predicted velocity at next time step (2 x Nx x Ny array)
% p = pressure (Nx x Ny array)
% h = dx and dy (2x1 array)
% dt = timestep
% ng = number of guardcells
% BC = boundary conditions (4x1 array)
% N = number of cells in each direction (2x1 array)
%
% returns:
% newP = pressure at next time step
	newP = p;
	div = divergence(uStar,vStar,h);
	rhs = div;
	newP = poissonSolve(rhs,p,h(1),BC,ng);
	
	% check that the pressure solve worked
	%[uStar,newP] = fillBC(uStar,newP,nDim,ng,N,BC);
	%remainder = checkPressureSolve(newP,rhs,h,op_choice(2),BC,ng,nDim)
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function [uNext,vNext] = correct(uStar,vStar,u,v,newP,h,dt,ng)
% correct the velocity with the new pressure field 
% this forces the velocity field to be divergence free
% see section XX of the documentation for more details
%
% inputs:
% uStar = predicted velocity at next time step (2 x Nx x Ny array)
% u = velocity (2 x Nx x Ny array)
% newP = pressure at next time step
% h = dx and dy (2x1 array)
% dt = timestep
% ng = number of guardcells
%
% returns:
% uNext = velocity at next time step
	dpdx = firstDeriv(newP,1,h,'pU');
	dpdy = firstDeriv(newP,2,h,'pV');
	
	uNext = uStar - dpdx*dt;
	vNext = vStar - dpdy*dt;
endfunction
