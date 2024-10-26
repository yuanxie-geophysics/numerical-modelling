%% This code is to calculate a density anormaly in mantle lithosphere
% Yuan Xie, 20.10.2022
% Homework 4

% 1) Clear memory and figures
clear
clf
close all
dbstop if error

% 2) Define Numerical model
xsize=100*1000; % Horizontal model size, m
ysize=100*1000; % Vertical model size, m
Nx=35; % Horizontal grid resolution
Ny=45; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m
x=0:dx:xsize; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize; % Vertical coordinates of basic grid points, m
xVx=0:dx:xsize; % Horizontal coordinates of vx grid points, m
yVx=-dy/2:dy:ysize-dy/2; % Vertical coordinates of vx grid points, m
xVy=-dx/2:dx:xsize-dx/2; % Horizontal coordinates of vy grid points, m
yVy=0:dy:ysize; % Vertical coordinates of vy grid points, m
xP=-dx/2:dx:xsize-dx/2; % Horizontal coordinates of P grid points, m
yP=-dy/2:dy:ysize-dy/2; % Vertical coordinates of P grid points, m
Vx=zeros(Ny1,Nx1); % Vx, m/s
Vy=zeros(Ny1,Nx1); % Vy, m/s
P=zeros(Ny1,Nx1); % Pressure, Pa
gy=10; % Gravity acceleration, m/s^2
rho=zeros(Ny,Nx); % Density, kg/m^3
eta=1e19; % Viscosity, Pa*s

% Define density field, rho
C=[xsize/2,ysize/2];% center point
radium=20*1000; % radium of density anormaly,m
for i=1:1:Ny
    for j=1:1:Nx
        px=x(j);
        py=y(i); % location of each point
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d <= radium
            rho(i,j)=3200;  % density
        else
            rho(i,j)=3300;  % density
        end
    end
end

% 3) Define global matrixes L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts
            
% 4) Composing global matrixes L(), R()
% Going through all points of the grid (2D meant two loops)
% composing respective equations

% First loop - horizontal index j
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index gVx, gVy, gP
        gVx=((j-1)*Ny1+i-1)*3+1; % Vx
        gVy=gVx+1; % Vy
        gP =gVx+2; % P
        
        % 4a) Eq. for Vx (ghost BC x-stokes)
        if(i==Ny1 || j==Nx1) % ghost nodes
            L(gVx,gVx)=1; % Left part
            R(gVx,1)=0; % Right part
        else
            if (i==1 || j==1 || j==Nx || i==Ny )                % Boundary Condition    
                L(gVx,gVx)=1; % Left part
                R(gVx,1)=0; % Right part
            else
                % Left part
                L(gVx,gVx-Ny1*3)=eta/dx^2; 
                L(gVx,gVx-3)=eta/dy^2; 
                L(gVx,gVx)=-2*eta/dx^2-2*eta/dy^2;
                L(gVx,gVx+3)=eta/dy^2; 
                L(gVx,gVx+Ny1*3)=eta/dx^2; 
                L(gVx,gP+Ny1*3)=-1/dx; % P2
                L(gVx,gP)=1/dx; % P1
                % Right part
                R(gVx,1)=0;
            end
        end
        
        % 4b) Eq. for Vy (ghost BC x-stokes)
        if( j==Nx1 || i==Ny1) % ghost nodes
            L(gVy,gVy)=1; % Left part
            R(gVy,1)=0; % Right part            
        else
            if (j==1 || i==1 || i==Ny  || j==Nx)    % Boundary Condition
                L(gVy,gVy)=1; % Left part
                R(gVy,1)=0; % Right part
            else
                % Left part
                L(gVy,gVy-Ny1*3)=eta/dx^2; 
                L(gVy,gVy-3)=eta/dy^2;
                L(gVy,gVy)=-2*eta/dy^2-2*eta/dx^2;
                L(gVy,gVy+3)=eta/dy^2; 
                L(gVy,gVy+Ny1*3)=eta/dx^2; 
                L(gVy,gP)=1/dy;  %P1
                L(gVy,gP+3)=-1/dy; %P2
                % Right part
                R(gVy,1)=-(rho(i,j-1)+rho(i,j))/2*gy;
            end
        end
        
        % 4c) Eq. for P (ghost BC Continuity-eq)
        if ( i==Ny1 || j==Nx1 )  % ghost nodes
            L(gP,gP)=1; % Left part
            R(gP,1)=0; % Right part
        else
            if ( i==1 || j==1 )   % Boundary Condition
                L(gP,gP)=1; % Left part
                R(gP,1)=0; % Right part
            else
                L(gP,gVx-Ny1*3)=-1/dx; %Vx1
                L(gP,gVx)=1/dx; % Vx2
                L(gP,gVy-3)=-1/dy; %Vy1
                L(gP,gVy)=1/dy; %Vy2
                R(gP,1)=0; % Right part
            end
        end
    end
end

% 5) Solving matrixes
S=L\R; 

% 6) Reload S--> PSI
for j=1:1:Nx1
    % Second loop - vertical index i
    for i=1:1:Ny1
        % Define global index g
        gVx=((j-1)*Ny1+i-1)*3+1; % Vx
        gVy=gVx+1; % Vy
        gP=gVx+2; % P
        % Reload solution
        Vx(i,j)=S(gVx);
        Vy(i,j)=S(gVy);
        P(i,j)=S(gP);
    end
end

% 7) Plot numerical solution

% colormap 
figure(1);colormap('Jet');
subplot(2,2,1)
pcolor(x,y,rho);% 
shading flat;
axis ij image;
colorbar
title('Density(kg/m^3)')

subplot(2,2,2)
pcolor(xP(1:Nx),yP(1:Ny),P(1:Ny,1:Nx))
shading interp;
axis ij image;
colorbar
title('Pressure(Pa)')

subplot(2,2,3)
pcolor(x,y,Vx(1:Ny,1:Nx))
shading interp;
axis ij image;
colorbar
title('vx(m/s)')

subplot(2,2,4)
pcolor(x,y,Vy(1:Ny,1:Nx))
shading interp;
axis ij image;
colorbar
title('vy(m/s)')


