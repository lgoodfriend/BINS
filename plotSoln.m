function plotSoln(u,v,ng,L,N,h)
% plot solution to BINS

	% u velocity--------------------------------------------------------------------------------------
	figure(1)
	xx = linspace(0,L,N+1);
	yy = linspace(h,L-h,N);
	[x,y] = meshgrid(xx,yy);
	surf(x,y,transpose(u(ng+1:ng+N+1,ng+1:ng+N)),'edgecolor','none')
	%caxis([-1 1])
	colorbar
	title('u velocity')
	view(2)
	grid off
	xlim([0 L]); ylim([0 L]);
	
	% v velocity--------------------------------------------------------------------------------------
	figure(2)
	xx = linspace(h,L-h,N);
	yy = linspace(0,L,N+1);
	[x,y] = meshgrid(xx,yy);
	surf(x,y,transpose(v(ng+1:ng+N,ng+1:ng+N+1 )),'edgecolor','none')
	%caxis([-1 1])
	colorbar
	title('v velocity')
	view(2)	
	grid off
	xlim([0 L]); ylim([0 L]);
	 
endfunction
