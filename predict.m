function [uStar,vStar] = predict(u,v,u_advection,v_advection,u_diffusion,v_diffusion,dt)
% estimates the velocity at the next time step
% the prediction step is discussed in Section 3.2, starting on page 5
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
end

