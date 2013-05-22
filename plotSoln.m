function plotSoln(u,v,ng,L,N,h)
	% plot solution to BINS
	% 
	% inputs:
	% u: x-direction velocity, size [N+2*ng+1   N+2*ng]
	% v: y-direction velocity, size [N+2*ng   N+2*ng+1]
	% ng: number of ghost cells
	% L: physical dimension in x and y 
	% N: number of interior points in x and y direction
	% h: grid spacing dx and dy
	%
	% returns:
	% nothing

	% limits of internal (non-ghostcell) data
        ll = ng+1; ul = ng+N;

	% u velocity--------------------------------------------------------------------------------------
	figure()
	% set up grid
	xx = linspace(0,L,N+1); % grid of x faces
	yy = linspace(h/2,L-h/2,N);
	[y,x] = meshgrid(yy,xx);
	% plot
	surf(x,y,u(ll:ul+1,ll:ul),'edgecolor','none')
	colorbar
	title('u velocity')
        xlabel('x'); ylabel('y');
	view(2)
	grid off
	xlim([0 L]); ylim([0 L]);
	
	% v velocity--------------------------------------------------------------------------------------
	figure()
	% set up grid
	xx = linspace(h/2,L-h/2,N); % grid of y faces
	yy = linspace(0,L,N+1);
	[y,x] = meshgrid(yy,xx);
        % plot
	surf(x,y,v(ll:ul,ll:ul+1),'edgecolor','none')
	colorbar
	title('v velocity')
        xlabel('x'); ylabel('y');
	view(2)	
	grid off
	xlim([0 L]); ylim([0 L]);
	
	% speed and streamlines---------------------------------------------------------------------------
	figure()
	% average velocities to cell centers
	midU = 0.5*(u(1:end-1,:      ) + u(2:end,:    ));
	midV = 0.5*(v(:      ,1:end-1) + v(:    ,2:end));
        % set up grid
	xx = linspace(h/2,L-h/2,N); % grid of cell centers
        yy = linspace(h/2,L-h/2,N);
        [y,x] = meshgrid(yy,xx);
        % plot contours of speed
        speed = sqrt(midU.^2 + midV.^2);
        contourf(yy,xx,transpose(speed(ll:ul,ll:ul)),'LineStyle','none')
        % set up streamlines
      	sx = linspace(h/2,L-h/2,N/4); % where streamlines start
      	sy = linspace(h/2,L-h/2,N/4);
      	[sxx,syy] = meshgrid(sx,sy);
	% plot streamlines
        hold on
      	h_fig = streamline(stream2(y,x,transpose(midU(ll:ul,ll:ul)),...
		transpose(midV(ll:ul,ll:ul)),syy,sxx,[0.1,300]));
        set(h_fig,'Color','white')
        hold off
        xlabel('x'); ylabel('y');
	title('speed with streamlines')
	colorbar
	 
	% vorticity---------------------------------------------------------------------------------------
        figure()
        % calculate vorticity
        dudy = (u(:,2:end)-u(:,1:end-1))/h; % du/dy at lower left cell corners
        dvdx = (v(2:end,:)-v(1:end-1,:))/h; % dv/dx at lower left cell corners
        vorticity = dvdx(ll-1:ul-1,ll:ul) - dudy(ll:ul,ll-1:ul-1); % just the inner (non-ghost) points
	% set up grid
        xx = linspace(0,L-h,N); % grid of lower left cell corners
        yy = linspace(0,L-h,N);
        [y,x] = meshgrid(yy,xx);
	% plot
	surf(x,y,vorticity,'edgecolor','none')
	colorbar
	title('vorticity')
        xlabel('x'); ylabel('y');
	view(2)	
	grid off
	xlim([0 L]); ylim([0 L]);
	
end














