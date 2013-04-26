function newP = pressure(uStar,vStar,p,h,dt,ng,BC,N,A)
% solves for the pressure using the requirement that the velocity
% field must be divergence free:
% laplacian(p) =  divergence(uStar) / dt
% the Poisson pressure solve is discussed in Section 3.2, starting on page 5
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
	newP = poissonSolve(rhs,h,BC,ng,N,A);
	
end
%--------------------------------------------------------------------------------------------------
% functions called by pressure:
% divergence
% poissonSolve (in own file, poissonSolve.m)
%
%--------------------------------------------------------------------------------------------------
function divOut = divergence(u,v,h)
	% calculates the divergence
	% for an overview of first derivatives, see Section 3.3.2 on page 7
	% for details, see Section 4.1 on page 10
	%
	% inputs:
	% u: x-direction velocity, size [N+2*ng+1   N+2*ng]
	% v: y-direction velocity, size [N+2*ng   N+2*ng+1]
	% h: spatial step
	%
	% returns:
	% divOut: divergence of u and v, size [N+2*ng   N+2*ng]
	
	dudx = (u(2:end,:)-u(1:end-1,:))/h;
	dvdy = (v(:,2:end)-v(:,1:end-1))/h;
	divOut = dudx + dvdy;
end

