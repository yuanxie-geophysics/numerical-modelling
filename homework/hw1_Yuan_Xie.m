%% Homework 1
% Yuan Xie, 29,09,2022

%% solving of 1D poission equation
% d2PHI/dx^2=2*x^2-x/2+exp(x);
% with finite difference
% on regular grid


%0) Clear variables and figures
clear % memory
clf % figure


% 1) Define the numerical model

xsize= 1; % Horizontal size of the model, m
Nx=101; % Number of grid points in the horizontal direction
dx=xsize/(Nx-1); % Grid step size,m
x=0:dx:xsize; % Horizontal cordinates at points, m

%2) Define global matrixes L(), R()
L=sparse(Nx,Nx); % left hand side coefficients(sparse matrix)
R=zeros(Nx,1);  % right hanside values

%3) Composing global matrixs L(), R()
% Going through all points of the grid
for j=1:1:Nx
	% define type of equation depending on index j
	if(j==1 || j==Nx)
		% BC equation : 1*PHT(j)=0
		% Left hand side 
		L(j,j)=1; % PHI(j) 
		% right hand side
		R(j,1)=0;
	else
		% poisson equation d2PHI/dx^2=1
		% stencil
		%     PHI(j-1)   PHI(j)    PHI(j+1)
		%------[j-1]------[j]-------[j+1]-----
		% Discretrized equation
		% Phi(j-1)-2*phi(j)+phi(j+1)/dx^2=1
		%
		%left hand side
		L(j,j-1)=1/dx^2;
		L(j,j)=-2/dx^2;
		L(j,j+1)=1/dx^2;
		% right hand side
		R(j,1)=2*x(j)^2-x(j)/2+exp(x(j));
	end
end

% 4) Solving of the equations
S=L\R;


%5) Reload solution to array phi
PHI=zeros(Nx,1);
for j=1:1:Nx
	PHI(j)=S(j);
end

%6) Visualize phi(i)
figure(1);clf
plot(x,PHI,'-o r')




% Comparing to analatical solution
Nxa=1001; % Number of points for analatical solution plot
dxa=xsize/(Nxa-1); % Grid step for the 
xa=0:dxa:xsize; % Coordinates of points for plotting the a..,m
% Analatical solution : PHI=x^2/2-x*L/2;
C1=1/xsize-(xsize^3)/6+(xsize^2)/12-exp(xsize)/xsize;
C2=-1;
PHIa=xa.^4/6-xa.^3/12+exp(xa)+C1.*xa+C2;
hold on
plot(xa,PHIa,'b')
