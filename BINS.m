% Basic Incompressible Navier-Stokes solver
% an educational tool
%
% written in Matlab, the prototyping language of engineers
% (actually written in Octave, the open source prototyping language of engineers) 
%
% Lauren Goodfriend, 2012

function BINS(N,L,dt,T,BC,IC_choice,nu,ng)
% solves the 2D incompressible Navier-Stokes equations
% for an overview of the algorithm used here, see Section 3.2 on page 4
%
% inputs: 
% N = number of cells in x and y 
% L = physical dimension in x and y 
% dt = timestep
% T = final time
% BC = boundary conditions (4x1 array)
% IC_choice = choice of initial conditions defined in IC.m
% nu = molecular viscosity
% ng = number of ghostcells 
%
% returns:
% nothing
	close all

	% establish initial conditions
	h = L./(N-1); 
	[u,v,p] = IC(L,N,ng,IC_choice);
	[u,v,p] = fillBC(u,v,p,ng,N,BC);

	t=0;
	while t < T
		% step forward in time

		% predict next step velocity with advection and diffusion terms (equation 1)
		[u_advection, v_advection] = advection(u,v,h,ng); % the advection term 
		[u_diffusion, v_diffusion] = diffusion(u,v,h,ng,BC,nu); % the diffusion term

		[uStar,vStar] = predict(u,v,u_advection,v_advection,u_diffusion,v_diffusion,dt);  % estimate next time step velocity
		[uStar,vStar,p] = fillBC(uStar,vStar,p,ng,N,BC); % enforce boundary conditions	

		% project the velocity field onto a divergence-free field to get the pressure (equation 2)
		newP = pressure(uStar,vStar,p,h,dt,ng,BC,N); % do the projection
		[uStar,vStar,newP] = fillBC(uStar,vStar,newP,ng,N,BC); % enforce boundary conditions

		% correct the velocity to be divergence-free using the updated pressure (equation 3)
		[newU,newV] = correct(uStar,vStar,newP,h,dt,ng);
		u=newU; v=newV;
		p=newP;
		[u,v,p] = fillBC(u,v,p,ng,N,BC);  % enforce boundary conditions

		% update the time
		t=t+dt;
	endwhile

	% plot the solution!
	plotSoln(u,v,ng,L,N,h)

                        % same the end time data in case we want it again
	save( ['./BINS_output',int2str(N(1)),'.mat'], ...
	             'u','v', 'p', 'ng', 'nu', 'L', 'N', 'T');

endfunction


