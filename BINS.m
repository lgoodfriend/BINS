% Berkeley Incompressible Navier-Stokes solver
% an educational tool
%
% written in Matlab, the prototyping language of engineers
% (actually written in Octave, the open source prototyping language of engineers) 
%
% Lauren Goodfriend, 2012

function BINS(N,L,dt,T,BC,IC_choice,nu,ng)
% solves the 2D incompressible Navier-Stokes equations
% for an overview of the algorithm used here, see Section 3.2 on page 5
%
% inputs: 
% N = number of cells in x and y 
% L = physical dimension in x and y 
% dt = timestep
% T = final time
% BC = boundary conditions, size [4] 
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
        % define Poisson operator for pressure solution
	A = make_matrix(N,ng,h,BC);

	t=0;
	while t < T
		% step forward in time

		% predict next step velocity with advection and diffusion terms (equations 1 and 2)
		[u_advection, v_advection] = advection(u,v,h); % the advection term 
		[u_diffusion, v_diffusion] = diffusionB(u,v,h,nu); % the diffusion term

		[uStar,vStar] = predict(u,v,u_advection,v_advection,u_diffusion,v_diffusion,dt);  % estimate next time step velocity
		[uStar,vStar,p] = fillBC(uStar,vStar,p,ng,N,BC); % enforce boundary conditions	

		% project the velocity field onto a divergence-free field to get the pressure (equation 3)
		newP = pressure(uStar,vStar,p,h,dt,ng,BC,N,A); % do the projection		
		[uStar,vStar,newP] = fillBC(uStar,vStar,newP,ng,N,BC); % enforce boundary conditions
		
		% correct the velocity to be divergence-free using the updated pressure (equations 4 and 5)
		[newU,newV] = correct(uStar,vStar,newP,h,ng,N,dt);
		u=newU; v=newV;
		p=newP;
		[u,v,p] = fillBC(u,v,p,ng,N,BC);  % enforce boundary conditions

		% update the time
		t=t+dt;
	end

	% plot the solution!
	plotSoln(u,v,ng,L,N,h)

        % save the end time data in case we want it again
        T = t;
	save( ['./BINS_output',int2str(N),'.mat'], ...
	             'u','v', 'p', 'ng', 'nu', 'L', 'N', 'T');

end


