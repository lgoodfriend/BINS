function [uOut,vOut,pOut] = fillBC(u,v,p,ng,N,BC);
	% fills the ghost cells of u, v, and p to enforce boundary conditions
	% 
	% inputs:
	% u: x-direction velocity (N+2*ng+1 x N+2*ng array)
	% v: y-direction velocity (N+2*ng x N+2*ng+1 array)
	% p: pressure (N+2*ng x N+2*ng array)
	% ng: number of ghost cells
	% N: number of interior points in x and y direction
	% BC: boundary condition tags in x and y direction (4x1 array)
	%
	% returns:
	% uOut: x-direction velocity with updated ghost cells (N+2*ng+1 x N+2*ng array)
	% vOut: y-direction velocity with updated ghost cells (N+2*ng x N+2*ng+1 array)
	% wOut: z-direction velocity with updated ghost cells (N+2*ng x N+2*ng array)
	uOut = u; vOut = v; pOut = p;
	
	for gcell=1:ng % loop over ghost cells
		% lower x BC-------------------------------------------------------------------------------------
		if BC(1)==sqrt(-1) % periodic, filled with upper x inner data
			uOut(gcell,:) = u(N+gcell+1,:);
			vOut(gcell,:) = v(N+gcell,:);
			pOut(gcell,:) = p(N+gcell,:);
		else % wall moving at speed BC parallel to itself
			uOut(gcell,:) = 0; uOut(gcell+1,:) = 0; 
			vOut(gcell,:) = 2*BC(1) - v(ng+1,:);
			pOut(gcell,:) = p(ng+1,:);
		endif
		% upper x BC-------------------------------------------------------------------------------------
		if BC(2)==sqrt(-1) % periodic, filled with lower x inner data
			uOut(N+ng+gcell+1,:) = u(ng+gcell,:);
			vOut(N+ng+gcell,:) = v(ng+gcell,:);
			pOut(N+ng+gcell,:) = p(ng+gcell,:);
		else % wall moving at speed BC parallel to itself
			uOut(N+ng+gcell,:) = 0; uOut(N+ng,:) = 0; 
			vOut(N+ng+gcell,:) = 2*BC(2) - v(ng+N,:);
			pOut(N+ng+gcell,:) = p(ng+N,:);
		endif
		% lower y BC-------------------------------------------------------------------------------------
		if BC(3)==sqrt(-1) % periodic, filled with upper y inner data
			uOut(:,gcell) = u(:,N+gcell);
			vOut(:,gcell) = v(:,N+gcell+1);
			pOut(:,gcell) = p(:,N+gcell);
		else % wall moving at speed BC parallel to itself
			uOut(:,gcell) = 2*BC(3) - u(:,ng+1);
			vOut(:,gcell) = 0; vOut(:,gcell+1) = 0;
			pOut(:,gcell) = p(:,ng+1);
		endif
		% upper y BC-------------------------------------------------------------------------------------
		if BC(4)==sqrt(-1) % periodic, filled with lower y inner data
			uOut(:,N+ng+gcell) = u(:,ng+gcell);
			vOut(:,N+ng+gcell+1) = v(:,ng+gcell);
			pOut(:,N+ng+gcell) = p(:,ng+gcell);
		else % wall moving at speed BC parallel to itself
			uOut(:,N+ng+gcell) = 2*BC(4) - u(:,N+gcell);
			vOut(:,N+ng+gcell) = 0; vOut(:,N+ng) = 0;
			pOut(:,N+ng+gcell) = p(:,N+gcell);
		endif
	endfor % loop over ghost cells

endfunction