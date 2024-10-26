%% This code is to calculate a density anormaly in mantle lithosphere
% With variable viscosity
% staggered nodes
% two phase flow
% no slip: = 1; free slip: = -1;
% Yuan Xie, 28.11.2022
% Homework 10

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
VxD=zeros(Ny1,Nx1); % Vx, m/s
VyD=zeros(Ny1,Nx1); % Vy, m/s
Pf=zeros(Ny1,Nx1); % Pressure, Pa
etaf=10;
rhof=2500;
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
%% Define fluid viscosity field, etafiPf
C=[xsize/2,ysize/2];% center point
etafiPf=zeros(Ny1,Nx1);
etafiM=10^21;
etafiP=10^20;
for i=1:1:Ny1
    for j=1:1:Nx1
        % compose ETAP on P nodes
        px=xP(j);
        py=yP(i);
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d <= radius
            etafiPf(i,j)=etafiP;
        else
            etafiPf(i,j)=etafiM;
        end
    end
end

%% Define fluid viscosity field, etafiPf
C=[xsize/2,ysize/2];% center point
kX=zeros(Ny1,Nx1);
kY=zeros(Ny1,Nx1);
kM=10^(-12);
kP=10^(-11);
for i=1:1:Ny1
    for j=1:1:Nx1
        % compose kfi on VxD nodes
        px=xVx(j);
        py=yVx(i); % location of each point
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d <= radius
            kX(i,j)=kP;
        else
            kX(i,j)=kM;
        end
        % compose kfi on VyD nodes
        px=xVy(j);
        py=yVy(i);
        d=sqrt((px-C(1))^2+(py-C(2))^2); % distance between center point and each point
        if d <= radius
            kY(i,j)=kP;
        else
            kY(i,j)=kM;
        end
    end
end
% 3) Define global matrixes L(), R()
N=Nx1*Ny1*6; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts

% 4) Composing global matrixes L(), R()
% Going through all points of the grid (2D meant two loops)
% composing respective equations

% Boundary conditions: free slip=-1; No Slip=1
bc=-1;

% First loop - horizontal index j
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index gVx, gVy, gP
        gVx =((j-1)*Ny1+(i-1))*6+1; % Vx,
        gVy =gVx+1; % Vy
        gP  =gVx+2; % P
        gVxD=gVx+3; % Vx
        gVyD=gVx+4;
        gPf =gVx+5;
        
        % 4a) Eq. for Vx (ghost BC x-stokes)
        if j==Nx1 % ghost nodes
            L(gVx,gVx)=1; % Left part
            R(gVx,1)=0; % Right part
        else
            if (i==1 || j==1 || j==Nx || i==Ny1 )    % Boundary Condition
                L(gVx,gVx)=1; % Left part
                R(gVx,1)=0; % Right part
                if (i==1)
                    L(gVx,gVx)=bc; % Left part
                    L(gVx,gVx+6)=1;
                end
                % Bottom boundary
                if (i==Ny1)
                    L(gVx,gVx-6)=bc; % Left part
                    L(gVx,gVx)=1;
                end
            else
                % viscosity
                eta1=eta(i-1,j);
                eta2=eta(i,j);
                etaP1=etaP(i,j);
                etaP2=etaP(i,j+1);
                % Left part
                L(gVx,gVx-Ny1*6)=2*etaP1/dx^2;  % vx1
                L(gVx,gVx-6)=eta1/dy^2;  % vx2
                L(gVx,gVx)=-2*(etaP1+etaP2)/dx^2-(eta1+eta2)/dy^2; % vx3
                L(gVx,gVx+6)=eta2/dy^2; % vx4
                L(gVx,gVx+Ny1*6)=2*etaP2/dx^2; % vx5
                L(gVx,gVy-6)=eta1/dx/dy; % vy1
                L(gVx,gVy)=-eta2/dx/dy; % vy2
                L(gVx,gVy+Ny1*6-6)=-eta1/dx/dy; % vy3
                L(gVx,gVy+Ny1*6)=eta2/dx/dy; % vy4
                L(gVx,gP+Ny1*6)=-1/dx; % P2
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
                    L(gVy,gVy+6*Ny1)=bc; % Left part
                    L(gVy,gVy)=1;
                end
                % Right boundary
                if (j==Nx1)
                    L(gVy,gVy-6*Ny1)=bc; % Left part
                    L(gVy,gVy)=1;
                end
                
            else
                % viscosity
                eta1=eta(i,j-1);
                eta2=eta(i,j);
                etaP1=etaP(i,j);
                etaP2=etaP(i+1,j);
                % Left part
                L(gVy,gVy-Ny1*6)=eta1/dx^2; % vy1
                L(gVy,gVy-6)=2*etaP1/dy^2; % vy2
                L(gVy,gVy)=-2*(etaP2+etaP1)/dy^2-(eta2+eta1)/dx^2; % vy3
                L(gVy,gVy+6)=2*etaP2/dy^2;  % vy4
                L(gVy,gVy+Ny1*6)=eta2/dx^2;  % vy5
                L(gVy,gVx-Ny1*6)=eta1/dx/dy; % vx1
                L(gVy,gVx-Ny1*6+6)=-eta1/dx/dy; % vx2
                L(gVy,gVx)=-eta2/dx/dy;  % vx3
                L(gVy,gVx+6)=eta2/dx/dy; % vx4
                L(gVy,gP)=1/dy;  %P1
                L(gVy,gP+6)=-1/dy; %P2
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
                L(gP,gVx-Ny1*6)=-1/dx; %Vx1
                L(gP,gVx)=1/dx; % Vx2
                L(gP,gVy-6)=-1/dy; %Vy1
                L(gP,gVy)=1/dy; %Vy2
                L(gP,gP)=1/etafiPf(i,j); %P
                L(gP,gPf)=-1/etafiPf(i,j); %Pf
                R(gP,1)=0; % Right part
            end
        end
        % 4d) Eq. for VxD
        if j==Nx1 % ghost nodes
            L(gVxD,gVxD)=1; % Left part
            R(gVxD,1)=0; % Right part
        else
            if (i==1 || j==1 || j==Nx || i==Ny1 )    % Boundary Condition
                L(gVxD,gVxD)=1; % Left part
                R(gVxD,1)=0; % Right part
                if (i==1)
                    L(gVxD,gVxD)=bc; % Left part
                    L(gVxD,gVxD+6)=1;
                end
                % Bottom boundary
                if (i==Ny1)
                    L(gVxD,gVxD-6)=bc; % Left part
                    L(gVxD,gVxD)=1;
                end
            else
                % Left part
                L(gVxD,gVxD)=etaf/kX(i,j);
                L(gVxD,gPf)=-1/dx; %Pf1
                L(gVxD,gPf+Ny1*6)=1/dx; % Pf2
                % Right part
                R(gVxD,1)=0;
            end
            % 4e) Eq. for VyD
            if i==Ny1 % ghost nodes
                L(gVyD,gVyD)=1; % Left part
                R(gVyD,1)=0; % Right part
            else
                if (j==1 || i==1 || i==Ny  || j==Nx1)    % Boundary Condition
                    L(gVyD,gVyD)=1; % Left part
                    R(gVyD,1)=0; % Right part
                    
                    % Left boundary
                    if (j==1)
                        L(gVyD,gVyD+6*Ny1)=bc; % Left part
                        L(gVyD,gVyD)=1;
                    end
                    % Right boundary
                    if (j==Nx1)
                        L(gVyD,gVyD-6*Ny1)=bc; % Left part
                        L(gVyD,gVyD)=1;
                    end
                else
                    
                    % left part
                    L(gVyD,gVyD)=etaf/kY(i,j); % VyD
                    L(gVyD,gPf)=-1/dy;  %Pf1
                    L(gVyD,gPf+6)=1/dy; %Pf2
                    % right part
                    R(gVyD,1)=rhof*gy;
                end
            end
            % 4f) Eq. for Pf
            if ( i==1 || j==1 || i==Ny1 || j==Nx1 )  % ghost nodes
                L(gPf,gPf)=1; % Left part
                R(gPf,1)=0; % Right part
            else
                L(gPf,gVxD-Ny1*6)=-1/dx; %VxD1
                L(gPf,gVxD)=1/dx; % VxD2
                L(gPf,gVyD-6)=-1/dy; %VyD1
                L(gPf,gVyD)=1/dy; %VyD2
                L(gPf,gP)=-1/etafiPf(i,j); %P
                L(gPf,gPf)=1/etafiPf(i,j); %Pf
                R(gPf,1)=0; % Right part
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
        gVx =((j-1)*Ny1+(i-1))*6+1; % Vx,
        gVy =gVx+1; % Vy
        gP  =gVx+2; % P
        gVxD=gVx+3; % Vx
        gVyD=gVx+4;
        gPf =gVx+5;
        
        % Reload solution
        Vx(i,j)=S(gVx);
        Vy(i,j)=S(gVy);
        P(i,j)=S(gP);
        VxD(i,j)=S(gVxD);
        VyD(i,j)=S(gVyD);
        Pf(i,j)=S(gPf);
    end
end

% 7) Plot numerical solution

% colormap
figure(1);colormap('Jet');
subplot(3,4,1)
pcolor(xVy,yVy,rho);%
shading flat;
axis ij image;
colorbar
title('Density(kg/m^3)')

subplot(3,4,2)
pcolor(x,y,log10(eta))
shading flat;
axis ij image;
colorbar
title('log ETAB(Pa.s)')

subplot(3,4,3)
pcolor(xP,yP,log10(etaP))
shading flat;
axis ij image;
colorbar
title('log ETAP(Pa.s)')

subplot(3,4,4)
pcolor(xP,yP,log10(etafiPf))
shading flat;
axis ij image;
colorbar
title('log ETAP(Pa.s)')

subplot(3,4,5)
pcolor(xP,yP,P)
shading interp;
axis ij image;
colorbar
title('Pressure(Pa)')

subplot(3,4,6)
pcolor(xVx,yVx,Vx)
shading interp;
axis ij image;
colorbar
title('vx(m/s)')

subplot(3,4,7)
pcolor(xVy,yVy,Vy)
shading interp;
axis ij image;
colorbar
title('vy(m/s)')

subplot(3,4,8)
pcolor(xVy,yVy,kY)
shading flat;
axis ij image;
colorbar
title('Vy-Permeability(m^2)')

subplot(3,4,9)
pcolor(xP,yP,Pf)
shading interp;
axis ij image;
colorbar
title('Pressure-Melts(Pa)')

subplot(3,4,10)
pcolor(xVx,yVx,VxD)
shading interp;
axis ij image;
colorbar
title('vxD(m/s)')

subplot(3,4,11)
pcolor(xVy,yVy,VyD)
shading interp;
axis ij image;
colorbar
title('vyD(m/s)')

subplot(3,4,12)
dp=P-Pf;
pcolor(xP,yP,dp)
shading interp;
axis ij image;
colorbar
title('Pressure-F(Pa)')

% FREE SLIP
aaa(1,1)=P(7,5);    %  1373979305.37810
aaa(2,1)=Vx(7,5);    % -1.83205164223993e-09
aaa(3,1)=Vy(7,5);    %  2.47614639386926e-09
aaa(4,1)=Pf(7,5);    %  1397981024.77070
aaa(5,1)=VxD(7,5);   % -1.24473387981050e-11
aaa(6,1)=VyD(7,5);   % -5.91274927480662e-10
format long e
disp(aaa)
