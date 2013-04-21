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

	% u velocity--------------------------------------------------------------------------------------
	figure(1)
	xx = linspace(0,L,N+1);
	yy = linspace(h/2,L-h/2,N);
	[y,x] = meshgrid(yy,xx);
	surf(x,y,u(ng+1:ng+N+1,ng+1:ng+N),'edgecolor','none')
	colorbar
	title('u velocity')
	view(2)
	grid off
	xlim([0 L]); ylim([0 L]);
	
	% v velocity--------------------------------------------------------------------------------------
	figure(2)
	xx = linspace(h/2,L-h/2,N);
	yy = linspace(0,L,N+1);
	[y,x] = meshgrid(yy,xx);
	surf(x,y,v(ng+1:ng+N,ng+1:ng+N+1),'edgecolor','none')
	colorbar
	title('v velocity')
	view(2)	
	grid off
	xlim([0 L]); ylim([0 L]);
	
	% speed and streamlines---------------------------------------------------------------------------
	figure(3)
	midU = 0.5*(u(1:end-1,:      ) + u(2:end,:    )); % velocities at cell centers
	midV = 0.5*(v(:      ,1:end-1) + v(:    ,2:end));
	xx = linspace(h/2,L-h/2,N); % grid
        yy = linspace(h/2,L-h/2,N);
        [y,x] = meshgrid(yy,xx);
        % contours of speed
        speed = sqrt(midU.^2 + midV.^2);
        contourf(yy,xx,transpose(speed(ng+1:ng+N,ng+1:ng+N)),'LineStyle','none')
        % streamlines
      	sx = linspace(h/2,L-h/2,N/4); % where streamlines start
      	sy = linspace(h/2,L-h/2,N/4);
      	[sxx,syy] = meshgrid(sx,sy);
        hold on
      	h = streamline(stream2(y,x,transpose(midU(ng+1:ng+N,ng+1:ng+N)),transpose(midV(ng+1:ng+N,ng+1:ng+N)),syy,sxx,[0.1,100]));
        set(h,'Color','white')
        hold off
	 
end
