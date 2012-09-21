function plotSoln(u,v,p,ng,L,N)
% plot solution to BINS
	% make x,y grids
	xx = linspace(0,L(1),N(1));
	yy = linspace(0,L(2),N(2));
	[x,y] = meshgrid(xx,yy);
	
	% u velocity--------------------------------------------------------------------------------------
	figure(1)
	surf(y,x,u(ng+1:ng+N(1),ng+1:ng+N(2) ),'edgecolor','none')
	%caxis([-1 1])
	colorbar
	title('u velocity')
	view(2)
	grid off
	xlim([0 L(1)]); ylim([0 L(2)]);
	
	% v velocity--------------------------------------------------------------------------------------
	figure(2)
	surf(y,x,v(ng+1:ng+N(1),ng+1:ng+N(2) ),'edgecolor','none')
	%caxis([-1 1])
	colorbar
	title('v velocity')
	view(2)	
	grid off
	xlim([0 L(1)]); ylim([0 L(2)]);
endfunction
