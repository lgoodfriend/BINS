function p = multigrid(f,pOld,h,N,BC)
        % This is unreliable and should not be used.  Sorry!


	% solves laplacian(p) = f
	% inputs:
	% f: right hand side of Poisson equation (N x N array)
	% h: grid spacing dx and dy
	% N: number of grid points in x and y
	% BC: array defining boundary conditions (4 x 1 array)
	% 
	% returns:
	% p: the solution to the Poisson equation, the new pressure (N x N array)
	p = pOld;
	eps = 1e-8; error = 1;
	iter = 0;
	not_done = true;
        errorTime = zeros(100,1);
	while not_done % loop through iterations of the multigrid cycle
		% one iteration
		p = mg(p,f,h,N,BC);
		iter = iter+1;
		% are we done?
		LaP = laplacian(p,h,N,BC);
		error = max(max(abs( LaP -f )));
		% if we're done, return
                errorTime(iter)=error;
		if error<=eps; return; end;
		if iter>100;
                   figure(); imagesc(LaP-f); colorbar; set(gca,'YDir','normal'); 
                   %figure(); imagesc(p); colorbar; set(gca,'YDir','normal'); 
                   %figure(); imagesc(LaP); colorbar; set(gca,'YDir','normal'); 
                   %figure(); imagesc(f); colorbar; set(gca,'YDir','normal'); 
                   figure(); semilogy(errorTime);
                   a-5/0; return; end;
	end
end
%--------------------------------------------------------------------------------
	%  multigrid calls:
	%	mg
	%	laplacian
	%	gauss_seidel
	%	average
	%	interpolate
        %       BC_mg
%--------------------------------------------------------------------------------
function p = mg(p,f,h,N,BC)
	% recursive function that goes up and down the multigrid "V"
	% inputs:
	% p: the current solution to the Poisson equation (N x N array)
	% f: the current right hand side to the Poisson equation (N x N array)
	% h: the current grid spacing dx and dy
	% N: the current number of grid points in x and y
	% BC: the array of boundary conditions (4 x 1 array)
	% 
	% returns:
	% p: updated solution to the Poisson equation

	% relax solution with Gauss-Seidel
	p = gauss_seidel(p,f,h,N,BC);
	p = gauss_seidel(p,f,h,N,BC);
	if N>1
		% calculate remainder
		LaP = laplacian(p,h,N,BC);
		r = LaP - f;
		% go down the "V"
		rc = average(r,N);
		delta = zeros(size(rc));
		rc = mg(delta,rc,2*h,N/2,BC);
		% go up the "V"
		r = interpolate(delta,N);
		p = p+r;
	end
	p = gauss_seidel(p,f,h,N,BC);
	p = gauss_seidel(p,f,h,N,BC);	
end
%--------------------------------------------------------------------------------
function LaP = laplacian(p,h,N,BC)
	% make domain with ghost cells
	pGC = zeros(N+2,N+2);
	pGC(2:N+1,2:N+1) = p;
	% fill domain's ghost cells
	pGC = BC_mg(pGC,1,N,BC);
	
        %pGC(1,:) = 0; pGC(:,1) = 0; pGC(end,:) = 0; pGC(:,end)=0;
	% calculate laplacian on inner domain
	LaP = zeros(N,N);
	LaP = (pGC(1:N,2:N+1) - 2*pGC(2:N+1,2:N+1) + pGC(3:N+2,2:N+1))/(h^2) ...% d^2p/dx^2
	    + (pGC(2:N+1,1:N) - 2*pGC(2:N+1,2:N+1) + pGC(2:N+1,3:N+2))/(h^2);   % d^2p/dy^2
end
%--------------------------------------------------------------------------------
function p = gauss_seidel(p,f,h,N,BC)
	% make domain with ghost cells
	pGC = zeros(N+2,N+2);
	pGC(2:N+1,2:N+1) = p;
	% fill domain's ghost cells
	pGC = BC_mg(pGC,1,N,BC);
	
        s = 1.6;
	% relax solution
        for i=2:N+1; for j=2:N+1;
		pGC(i,j) = s*0.25*( ...
                   pGC(i-1,j) + pGC(i+1,j) + ... % left + right +
       	           pGC(i,j-1) + pGC(i,j+1) + ... % down + up +
               	  -(h^2)*f(i-1,j-1) ) + ...      % rhs
		   (1-s)*pGC(i,j);
        end; end
	p = pGC(2:end-1,2:end-1);
	
end
%--------------------------------------------------------------------------------
function rc = average(r,N)
	% average solutions in a 2x2 region to approximate the solution at the next
	% smaller grid
	rc = zeros(N/2,N/2);
	rc = 0.25*( r(2:2:N,2:2:N) + r(1:2:N-1,2:2:N) + r(2:2:N,1:2:N-1) + r(1:2:N-1,1:2:N-1) );
end
%--------------------------------------------------------------------------------
function r = interpolate(delta,N)
	% inject solutions from the coarse to the fine grid
	r = zeros(N,N);
	r(2:2:N  ,2:2:N  ) = delta(1:N/2,1:N/2);
	r(1:2:N-1,2:2:N  ) = delta(1:N/2,1:N/2);
	r(2:2:N  ,1:2:N-1) = delta(1:N/2,1:N/2);
	r(1:2:N-1,1:2:N-1) = delta(1:N/2,1:N/2);

        % spline interpolate from the coarse to the fine grid
        % this has no effect on the convergence
        %xxC = 1:N/2; yyC = 1:N/2; [xC,yC] = ndgrid(xxC,yyC); % (x,y) coordinates of coarse grid
        %xxF = 0.75:0.5:N/2+0.25; yyF = 0.75:0.5:N/2+0.25; [xF,yF] = ndgrid(xxF,yyF); % (x,y) coordinates of fine grid
        %F = griddedInterpolant(xC,yC,delta,'spline'); % create interpolation function for delta
        %r = F(xF,yF); % interpolate onto fine grid 
end
%--------------------------------------------------------------------------------
function pOut = BC_mg(p,ng,N,BC);
	pOut = p;
	
	for gcell=1:ng % loop over ghost cells
		% lower x BC-------------------------------------------------
		if BC(1)==sqrt(-1) % periodic, filled with upper x inner data
			pOut(gcell,:) = p(N+gcell,:);
		else % wall moving at speed BC parallel to itself
			pOut(gcell,:) = p(ng+1,:);
		end
		% upper x BC-------------------------------------------------
		if BC(2)==sqrt(-1) % periodic, filled with lower x inner data
			pOut(N+ng+gcell,:) = p(ng+gcell,:);
		else % wall moving at speed BC parallel to itself
			pOut(N+ng+gcell,:) = p(ng+N,:);
		end
		% lower y BC-------------------------------------------------
		if BC(3)==sqrt(-1) % periodic, filled with upper y inner data
			pOut(:,gcell) = p(:,N+gcell);
			pOut(gcell,gcell) = pOut(gcell+N,gcell+N);
		else % wall moving at speed BC parallel to itself
			pOut(:,gcell) = p(:,ng+1);
                        pOut(gcell,gcell) = p(ng+1,ng+1);
		end
		% upper y BC-------------------------------------------------
		if BC(4)==sqrt(-1) % periodic, filled with lower y inner data
			pOut(:,N+ng+gcell) = p(:,ng+gcell);
		else % wall moving at speed BC parallel to itself
			pOut(:,N+ng+gcell) = p(:,N+gcell);
		end
	end % loop over ghost cells
end
%--------------------------------------------------------------------------------
