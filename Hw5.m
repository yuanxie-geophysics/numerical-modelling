%% This code is to calculate a density anormaly in mantle lithosphere
% With variable viscosity
% staggered nodes
% no slip: = 1; free slip: = -1;  
% Yuan Xie, 24.10.2022
% Homework 5

% 1) Clear memory and figures
clear all
clf

% 2) Define Numerical model
xsize=100*1000; % Horizontal model size, m
ysize=100*1000; % Vertical model size, m
Nx=35; % Horizontal grid resolution
Ny=45; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m
x=0:dx:xsize+dx; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize+dy; % Vertical coordinates of basic grid points, m
xVx=0:dx:xsize+dx; % Horizontal coordinates of vx grid points, m
yVx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
xVy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yVy=0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
xP=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yP=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m
Vx=zeros(Ny1,Nx1); % Vx, m/s
Vy=zeros(Ny1,Nx1); % Vy, m/s
P=zeros(Ny1,Nx1); % Pressure, Pa
gy=10; % Gravity acceleration, m/s^2
rho=zeros(Ny1,Nx1); % Density, kg/m^3
eta=zeros(Ny1,Nx1); % Viscosity, Pa*s
etaP=zeros(Ny1,Nx1); % Viscosity, Pa*s

% Define density field, rho
C=[xsize/2,ysize/2];% center point
radius=20*1000; % radium of density anormaly,m
for i=1:1:Ny1
    for j=1:1:Nx1
        px=xVy(j);
        py=yVy(i); % location of each point
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d <= radius
            rho(i,j)=3200;  % density
        else
            rho(i,j)=3300;  % density
        end
    end
end

% Define viscosity field, eta etaP
C=[xsize/2,ysize/2];% center point
etaIN=1e18;
etaOU=1e19;
for i=1:1:Ny1
    for j=1:1:Nx1
        % compose ETA on basic nodes
        px=x(j);
        py=y(i); % location of each point
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d <= radius
            eta(i,j)=etaIN;
        else
            eta(i,j)=etaOU;
        end
        % compose ETAP on P nodes
        pxP=xP(j);
        pyP=yP(i);
        d=sqrt((pxP-C(1))^2+(pyP-C(2))^2); % distance between center point and each point
        if d <= radius
            etaP(i,j)=etaIN;
        else
            etaP(i,j)=etaOU;
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

% Boundary conditions: free slip=-1; No Slip=1
bc=1;


% First loop - horizontal index j
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index gVx, gVy, gP
        gVx=((j-1)*Ny1+(i-1))*3+1; % Vx
        gVy=gVx+1; % Vy
        gP =gVx+2; % P
        
        % 4a) Eq. for Vx (ghost BC x-stokes)
        if j==Nx1 % ghost nodes
            L(gVx,gVx)=1; % Left part
            R(gVx,1)=0; % Right part
        else
            if (i==1 || j==1 || j==Nx || i==Ny1 )    % Boundary Condition
                L(gVx,gVx)=1; % Left part
                R(gVx,1)=0; % Right part
                % Right boundary
                if (j==Nx)
                    L(gVx,gVy)=1;
                    L(gVx,gVy+Ny1*3)=bc;
                end
                % Left boundary
                if (j==1)
                    L(gVx,gVy)=1;
                    L(gVx,gVy+Ny1*3)=bc;
                end
                % Top boundary
                if (i==1)
                    L(gVx,gVx)=bc; % Left part
                    L(gVx,gVx+3)=1;
                end
                % Bottom boundary
                if (i==Ny1)
                    L(gVx,gVx-3)=bc; % Left part
                    L(gVx,gVx)=1;
                end
            else
                % viscosity
                eta1=eta(i-1,j);
                eta2=eta(i,j);
                etaP1=etaP(i,j);
                etaP2=etaP(i,j+1);
                % Left part
                L(gVx,gVx-Ny1*3)=2*etaP1/dx^2;  % vx1
                L(gVx,gVx-3)=eta1/dy^2;  % vx2
                L(gVx,gVx)=-2*(etaP1+etaP2)/dx^2-(eta1+eta2)/dy^2; % v3
                L(gVx,gVx+3)=eta2/dy^2; % vx4
                L(gVx,gVx+Ny1*3)=2*etaP2/dx^2; % vx5
                L(gVx,gVy-3)=eta1/dx/dy; % vy1
                L(gVx,gVy)=-eta2/dx/dy; % vy2
                L(gVx,gVy+Ny1*3-3)=-eta1/dx/dy; % vy3
                L(gVx,gVy+Ny1*3)=eta2/dx/dy; % vy4
                L(gVx,gP+Ny1*3)=-1/dx; % P2
                L(gVx,gP)=1/dx; % P1
                % Right part
                R(gVx,1)=0;
            end
        end
        % 4b) Eq. for Vy (ghost BC x-stokes)
        if i==Ny1 % ghost nodes
            L(gVy,gVy)=1; % Left part
            R(gVy,1)=0; % Right part
        else
            if (j==1 || i==1 || i==Ny  || j==Nx1)    % Boundary Condition
                L(gVy,gVy)=1; % Left part
                R(gVy,1)=0; % Right part
                % Left boundary
                if (j==1)
                    L(gVy,gVy+3*Ny1)=bc; % Left part
                    L(gVy,gVy)=1;
                end
                % Right boundary
                if (j==Nx1)
                    L(gVy,gVy-3*Ny1)=bc; % Left part
                    L(gVy,gVy)=1;
                end
                % Top boundary
                if (i==1)
                    L(gVy,gVx+3)=bc; % Left part
                    L(gVy,gVx)=1;
                end
                % Bottom boundary
                if (i==Ny)
                    L(gVy,gVx+3)=bc; % Left part
                    L(gVy,gVx)=1;
                end
            else
                % viscosity
                eta1=eta(i,j-1);
                eta2=eta(i,j);
                etaP1=etaP(i,j);
                etaP2=etaP(i+1,j);
                % Left part
                L(gVy,gVy-Ny1*3)=eta1/dx^2; % vy1
                L(gVy,gVy-3)=2*etaP1/dy^2; % vy2
                L(gVy,gVy)=-2*(etaP2+etaP1)/dy^2-(eta2+eta1)/dx^2; % vy3
                L(gVy,gVy+3)=2*etaP2/dy^2;  % vy4
                L(gVy,gVy+Ny1*3)=eta2/dx^2;  % vy5
                L(gVy,gVx-Ny1*3)=eta1/dx/dy; % vx1
                L(gVy,gVx-Ny1*3+3)=-eta1/dx/dy; % vx2
                L(gVy,gVx)=-eta2/dx/dy;  % vx3
                L(gVy,gVx+3)=eta2/dx/dy; % vx4
                L(gVy,gP)=1/dy;  %P1
                L(gVy,gP+3)=-1/dy; %P2
                % Right part
                R(gVy,1)=-rho(i,j)*gy;
            end
        end
        
        % 4c) Eq. for P (ghost BC Continuity-eq)
        if ( i==1|| j==1 || i==Ny1 || j==Nx1 )  % ghost nodes
            L(gP,gP)=1; % Left part
            R(gP,1)=0; % Right part
        else
            if ( i==2 && j==2 )   % Boundary Condition
                L(gP,gP)=1; % Left part
                R(gP,1)=1e9; % Right part
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
        gVx=((j-1)*Ny1+(i-1))*3+1; % Vx
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
subplot(2,3,1)
pcolor(xVy,yVy,rho);% 
shading flat;
axis ij image;
colorbar
title('Density(kg/m^3)')

subplot(2,3,2)
pcolor(x,y,log10(eta))
shading flat;
axis ij image;
colorbar
title('log ETAB(Pa.s)')

subplot(2,3,3)
pcolor(xP,yP,log10(etaP))
shading flat;
axis ij image;
colorbar
title('log ETAP(Pa.s)')

subplot(2,3,4)
pcolor(xP,yP,P)
shading interp;
axis ij image;
colorbar
title('Pressure(Pa)')

subplot(2,3,5)
pcolor(xVx,yVx,Vx)
shading interp;
axis ij image;
colorbar
title('vx(m/s)')

subplot(2,3,6)
pcolor(xVy,yVy,Vy)
shading interp;
axis ij image;
colorbar
title('vy(m/s)')


% FREE SLIP 
aaa(1,1)=P(7,5);  %  1374927109.05287
aaa(2,1)=Vx(7,5);    % -1.75951244594608e-09
aaa(3,1)=Vy(7,5);   %  1.79677698022029e-09 
% NO SLIP 
% aaa(1,1)=pr(7,5);    %  1374959400.70958
% aaa(2,1)=vx(7,5);    % -5.44992100189968e-10
% aaa(3,1)=vy(7,5);    %  5.87771380138805e-10
if (bc==1)
    disp('No slip');
else 
    disp('Free slip');
end
format long e
disp(aaa)