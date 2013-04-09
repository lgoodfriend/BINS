function [uNext,vNext] = correct(uStar,vStar,newP,h,dt)
% correct the velocity with the new pressure field 
% this forces the velocity field to be divergence free
% the correction step is described in Section 3.2, starting on page 5
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

end
%--------------------------------------------------------------------------------------------------
% functions called by correct:
% correct_dpdx
% correct_dpdy
%
% both of the pressure derivative functions
% like the advection terms, these are very similar but differ due to the staggered grid
% for an overview of the first derivatives, see Section 3.3.2 on page 7
%--------------------------------------------------------------------------------------------------
function dpdx = correct_dpdx(p,h)
	% calculates dpdx on cell faces from p on cell centers
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% p: pressure, size [N+2*ng   N+2*ng]
	% h: spatial step
	%
	% returns:
	% dpdx: derivative of pressure in the x direction, size [N+2*ng+1   N+2*ng]
	
	% initialize to correct size
	size_p = size(p);
	dpdx = zeros( size_p(1)+1, size_p(2));
	% take derivative
	dpdx(2:end-1,:) = (p(2:end,:)-p(1:end-1,:))/h;
end

function dpdy = correct_dpdy(p,h)
	% calculates dpdx on cell faces from p on cell centers
	% for details, see Section 4.1 (page 10)
	%
	% inputs:
	% p: pressure, size [N+2*ng   N+2*ng]
	% h: spatial step
	%
	% returns:
	% dpdy: derivative of pressure in the y direction, size  [N+2*ng  N+2*ng+1]
	
	% initialize to correct size
	size_p = size(p);
	dpdy = zeros( size_p(1), size_p(2)+1);
	% take derivative
	dpdy(:,2:end-1) = (p(:,2:end)-p(:,1:end-1))/h;
end
