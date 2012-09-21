% calculus
% this file contains routines that estimate derivatives:
%	firstDeriv
%	secondDeriv
%	divergence

1; % this is needed in octave for the interpreter to find all the functions
%------------------------------------------------------------------------------------------------------------------------------------------
function derivOut = firstDeriv(u,dir,h,use)
	% calculate first derivative for various input and output configurations
	%
	% inputs:
	% u: variable to take derivative of
	% dir: direction of derivative, 1 for x direction, 2 for y direction
	% h: spatial step
	% use: tag marking destination of output derivative
	%	This is used to mark the dimensions of the input and output data.
	%	This is needed to allow derivatives of face centered data to be
	%	output to cell centers, and vice-versa.
	%
	% returns:
	% derivOut: first derivative of input variable u
	
	% find dimensions of input data iX, iY
	sizeInput = size(u);
                    iX = sizeInput(1); iY = sizeInput(2); 
                    
                    % declare dimensions of output data oX, oY
                    % lower limits of derivative data lX,lY
                    % and upper limits of derivative data uX,uY
                    if or(strcmp(use,'advU'),  strcmp(use,'pU')) % input data is NxN, output data is N+1 x N
                    	oX = iX+1; oY = iY; 
                    	lX=2; lY=2;
                    	uX = iX; uY = iY;
                    elseif or(strcmp(use,'advV'), strcmp(use,'pV')) % input data is NxN, output data is N x N+1
                    	oX = iX; oY = iY+1;
                    	lX=2; lY=2;
                    	uX = iX; uY = iY;
                    elseif strcmp(use,'divU') % input data is N+1 x N, output data is NxN
                    	oX = iX-1; oY = iY;
                    	lX=1;lY=2;
                    	uX = oX; uY = oY;
                    elseif strcmp(use,'divV') % input data is N x N+1, output data is NXN
                    	oX = iX; oY = iY-1;
                    	lX=2; lY=1;
                    	uX = oX; uY = oY;
                    endif
                    derivOut = zeros(oX,oY);
                    
                    % take a derivative
                    if dir==1 % d/dx
                    	derivOut(lX:uX,1:uY) = (u(2:end,:)-u(1:end-1,:))/h(1);
                    elseif dir==2 % d/dy
                    	derivOut(1:uX,lY:uY) = (u(:,2:end)-u(:,1:end-1))/h(2);
                    endif
                   		
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function derivOut = secondDeriv(u,dir,h)
	% calculate second derivative
	% This is much less complicated than the firstDeriv function because all output derivatives are 
	% located in the same part of the cell as the input variable.
	%
	% inputs:
	% u: variable to take derivative of
	% dir: direction of derivative, 1 for x direction, 2 for y direction
	% h: spatial step 
	%
	% returns:
	% derivOut: second derivative of input variable u

	derivOut = u;
	switch dir
		case{1} % x direction
			derivOut(2:end-1,:) = (u(1:end-2,:) + u(3:end,:) - 2*u(2:end-1,:))/(h(1)^2);
		case{2} % y direction
			derivOut(:,2:end-1) = (u(:,1:end-2) + u(:,3:end) - 2*u(:,2:end-1))/(h(2)^2);
	endswitch % dir
endfunction

%------------------------------------------------------------------------------------------------------------------------------------------
function divOut = divergence(u,v,h)
	% calculates the divergence
	%
	% inputs:
	% u: x-direction velocity
	% v: y-direction velocity
	% h: spatial step
	%
	% returns:
	% divOut: divergence of u and v
	
	dudx = firstDeriv(u,1,h,'divU');
	dvdy = firstDeriv(v,2,h,'divV');
	divOut = dudx + dvdy;
endfunction
