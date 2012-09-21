function [uOut,vOut,pOut] = fillBC(u,v,p,ng,N,BC);
	% fills the ghost cells of u, v, and p to enforce boundary conditions
	% 
	% inputs:
	% u: x-direction velocity
	% v: y-direction velocity
	% p: pressure
	% ng: number of ghost cells
	% N: number of interior points in x and y direction (2x1 array)
	% BC: boundary condition tags in x and y direction (4x1 array)
	%
	% returns:
	% uOut: x-direction velocity with updated ghost cells
	% vOut: y-direction velocity with updated ghost cells
	% wOut: z-direction velocity with updated ghost cells
	uOut = u; vOut = v; pOut = p;
	
	for gcell=1:ng % loop over ghost cells
		% lower x BC-------------------------------------------------------------------------------------
		if BC(1)==sqrt(-1) % periodic, filled with upper x inner data
			uOut(gcell,:) = u(N(1)+gcell+1,:);
			vOut(gcell,:) = v(N(1)+gcell,:);
			pOut(gcell,:) = p(N(1)+gcell,:);
		else % wall moving at speed BC parallel to itself
			uOut(gcell,:) = 0; uOut(gcell+1,:) = 0; 
			vOut(gcell,:) = 2*BC(1) - v(ng+1,:);
			pOut(gcell,:) = p(ng+1,:);
		endif
		% upper x BC-------------------------------------------------------------------------------------
		if BC(2)==sqrt(-1) % periodic, filled with lower x inner data
			uOut(N(1)+ng+gcell+1,:) = u(ng+gcell,:);
			vOut(N(1)+ng+gcell,:) = v(ng+gcell,:);
			pOut(N(1)+ng+gcell,:) = p(ng+gcell,:);
		else % wall moving at speed BC parallel to itself
			uOut(N(1)+ng+gcell,:) = 0; uOut(N(1)+ng,:) = 0; 
			vOut(N(1)+ng+gcell,:) = 2*BC(2) - v(ng+N(1),:);
			pOut(N(1)+ng+gcell,:) = p(ng+N(1),:);
		endif
		% lower y BC-------------------------------------------------------------------------------------
		if BC(3)==sqrt(-1) % periodic, filled with upper y inner data
			uOut(:,gcell) = u(:,N(2)+gcell);
			vOut(:,gcell) = v(:,N(2)+gcell+1);
			pOut(:,gcell) = p(:,N(2)+gcell);
		else % wall moving at speed BC parallel to itself
			uOut(:,gcell) = 2*BC(3) - u(:,ng+1);
			vOut(:,gcell) = 0; vOut(:,gcell+1) = 0;
			pOut(:,gcell) = p(:,ng+1);
		endif
		% upper y BC-------------------------------------------------------------------------------------
		if BC(4)==sqrt(-1) % periodic, filled with lower y inner data
			uOut(:,N(2)+ng+gcell) = u(:,ng+gcell);
			vOut(:,N(2)+ng+gcell+1) = v(:,ng+gcell);
			pOut(:,N(2)+ng+gcell) = p(:,ng+gcell);
		else % wall moving at speed BC parallel to itself
			uOut(:,N(2)+ng+gcell) = 2*BC(4) - u(:,N(2)+gcell);
			vOut(:,N(2)+ng+gcell) = 0; vOut(:,N(2)+ng) = 0;
			pOut(:,N(2)+ng+gcell) = p(:,N(2)+gcell);
		endif
	endfor % loop over ghost cells

endfunction