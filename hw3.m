%% This code is to calculate a density anormaly in mantle lithosphere
% Yuan Xie, 11.10.2022
% Homework 3
%%

% with finite differences
% on a regular grid

% 1) Clear memory and figures
clear
clf

% 2) Define numerical model
xsize=100*1000; % Horizontal model size, m
ysize=100*1000; % Vertical model size, m
Nx=35; % Horizontal resolution
Ny=45; % Vertical resolution
dx=xsize/(Nx-1); % Horizontal grid step,m
dy=ysize/(Ny-1); % Vertical grid step, m
x=0:dx:xsize; % Horizontal coordinates of grid points, m
y=0:dy:ysize;  % Vertical coordinates of grid points, m
eta=10^19; % Viscosity
gy=10; % gravity accelaration,m/s^2


% 3) Define global matrixes L(), R()
N=Nx*Ny; % Number of unknowns
L=sparse(N,N); % Left hand side coefficents
R=zeros(N,1); % Right hand side values
rho=zeros(Ny,Nx);

% ------------ CALCULATING OMEGA(vorticity, 1/s) ---------------
% 4) Composing of global matrixes L(),R()
% Going through all points of the grid (2D meant two loops)
% define density field, rho
rhoi=3200; % inside density, kg/m^3
rhoo=3300; % outside density, kg/m^3
radium=20*1000; % radium of density anormaly,m
C=[xsize/2,ysize/2];% center point
for j=1:1:Nx
    for i=1:1:Ny
        px=dx*(j-1);
        py=dy*(i-1); % location of each point
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d < radium
            rho(i,j)=rhoi;
        else
            rho(i,j)=rhoo;
        end
    end
end

% First loop - horizontal index j
for j=1:1:Nx
	% Second loop - vertical index i
	for i=1:1:Ny
		% Define global index g
		g=(j-1)*Ny+i;
		% Decide which equation to solve BC or Poisson eq
		if(j==1 || j==Nx || i==1 || i==Ny)
			% BS: OMEGA(i,j)=0 => 1*S(g)=0 
			L(g,g)=1; % Left hand side for FI(i,j)
			R(g,1)=0; % Right hand side 
			else
			% Left hand side (5 unknowns)
			L(g,g-Ny)=1/dx^2; %OMEGA_A
			L(g,g-1)=1/dy^2; %OMEGA_B
			L(g,g)=-2/dx^2-2/dy^2; %OMEGA_C
			L(g,g+1)=1/dy^2; %OMEGA_D
			L(g,g+Ny)=1/dx^2; %OMEGA_E
			% Right hand side
			R(g,1)=1/eta*gy*((rho(i,j+1)-rho(i,j-1))/(2*dx));
		end
	end
end

% 5) Solve global matrixes
S=L\R;
% 6) Reload S--> OMEGA
OMEGA=zeros(Ny,Nx); % Creat 
for j=1:1:Nx
	% Second loop - vertical index i
	for i=1:1:Ny
		% Define global index g
		g=(j-1)*Ny+i;
		% Reload S(g)->PHI
		OMEGA(i,j)=S(g);
	end
end

% ------------ CALCULATING PSI(STREAM FUNCTION, m^2/s) ---------------
% 4) Composing of global matrixes L(),R()
% Going through all points of the grid (2D meant two loops)
% First loop - horizontal index j
for j=1:1:Nx
	% Second loop - vertical index i
	for i=1:1:Ny
		% Define global index g
		g=(j-1)*Ny+i;
		% Decide which equation to solve BC or Poisson eq
		if(j==1 || j==Nx || i==1 || i==Ny)
			% BS: PSI(i,j)=0 => 1*S(g)=0 
			L(g,g)=1; % Left hand side for FI(i,j)
			R(g,1)=0; % Right hand side 
			else
			% Left hand side (5 unknowns)
			L(g,g-Ny)=1/dx^2; %OMEGA_A
			L(g,g-1)=1/dy^2; %OMEGA_B
			L(g,g)=-2/dx^2-2/dy^2; %OMEGA_C
			L(g,g+1)=1/dy^2; %OMEGA_D
			L(g,g+Ny)=1/dx^2; %OMEGA_E
			% Right hand side
			R(g,1)=OMEGA(i,j);
		end
	end
end

% 5) Solve global matrixes
S=L\R;
% 6) Reload S--> PSI
PSI=zeros(Ny,Nx); % Creat 
for j=1:1:Nx
	% Second loop - vertical index i
	for i=1:1:Ny
		% Define global index g
		g=(j-1)*Ny+i;
		% Reload S(g)->PHI
		PSI(i,j)=S(g);
	end
end

% ------------ CALCULATING VELOCITY Vx,Vy(m/s) ---------------
vy=zeros(Ny,Nx);
vx=zeros(Ny,Nx);
for j=1:1:Nx
    for i=1:1:Ny
        if(j==1 || j==Nx || i==1 || i==Ny)
        %BC: vx=0;vy=0
            vx(i,j)=0;
            vy(i,j)=0;
        else      
            vx(i,j)=(PSI(i+1,j)-PSI(i-1,j))/(2*dy);
            vy(i,j)=-(PSI(i,j+1)-PSI(i,j-1))/(2*dx);
        end
    end
end



% 7) Plot numerical solution

% colormap 
figure(1); clf;colormap('Jet');
subplot(2,3,1)
colormap;
pcolor(x,y,rho)
colorbar
title('density,kg/m^3')
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down	

subplot(2,3,2)
colormap;
pcolor(x,y,OMEGA)
colorbar
shading interp
title('vorticity,1/s')
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down	

subplot(2,3,3)
colormap;
pcolor(x,y,PSI)
colorbar
shading interp
title('stream function,m^2/s')
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down	

subplot(2,3,4)
colormap;
pcolor(x,y,vx)
colorbar
shading interp
title('vx-velocity,m/s')
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down	

subplot(2,3,5)
colormap;
pcolor(x,y,vy)
colorbar
shading interp
title('vy-velocity,m/s')
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down	

subplot(2,3,6)
colormap;
quiver(x,y,vx,vy,'black')
title('Velocity field, m/s')
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down	