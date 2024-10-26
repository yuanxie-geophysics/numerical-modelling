%% Homework 2_1
% Yuan Xie, 06,10,2022

%% Solving of 2D poission equation using Gauss-Seidel Iteration
% d2PHI/dx^2+d2PHI/dy^2=1
% with finite difference
% on regular grid


% 0) Clear variables and figures
clf;clear;close all;

% 1) Define numerical model
xsize=1; % Horizontal model size, m
ysize=1; % Vertical model size, m
Nx=35; % Horizontal resolution
Ny=45; % Vertical resolution
dx=xsize/(Nx-1); % Horizontal grid step,m
dy=ysize/(Ny-1); % Vertical grid step, m
x=0:dx:xsize; % Horizontal coordinates of grid points, m
y=0:dy:ysize;  % Vertical coordinates of grid points, m

% 2) Define global matrixes
PHI0=zeros(Ny,Nx);
PHI1=zeros(Ny,Nx);
R=zeros(Ny,Nx);
for i=1:1:Ny
  for j=1:1:Nx
    %Boundary nodes
    if(i==1 || i==Ny || j==1 || j==Nx)
        R(i,j)=0;
    else
    %Internal nodes
        R(i,j)=1;
    end
  end
end

% 3) Iteration loop
delta=1.5; % relaxation fator
Nmax=10;
for niter=1:1:Nmax
    for i=1:1:Ny
        for j=1:1:Nx
            if(j==1 || j==Nx || i==1 || i==Ny)
                PHI1(i,j)=0;
            else
                dR=R(i,j)-((PHI0(i,j-1)-2*PHI0(i,j)+PHI0(i,j+1))/(dx^2)+(PHI0(i-1,j)-2*PHI0(i,j)+PHI0(i+1,j))/(dy^2));
                PHI1(i,j)=PHI0(i,j)+dR/(-2/dx^2-2/dy^2)*delta;                
            end
            PHI0(i,j)=PHI1(i,j);
        end
    end
end

% 4) Reload PHI
PHI=PHI1;

% 5) Visualization
% colormap
figure(1); clf;
% Colormap
subplot(1,2,1)
pcolor(x,y,PHI)
colorbar
shading interp
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down
% 3D surface
subplot (1,2,2)
surf(x,y,-PHI)
colorbar
light
lighting phong
shading interp

