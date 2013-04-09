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
	 
end
