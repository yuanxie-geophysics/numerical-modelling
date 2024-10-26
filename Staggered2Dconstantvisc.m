% Solving of momentum and continuity equations
% with FD on staggered grid

% Clear variables
clear all

% Define Numerical model
xsize=100000; % Hortizontal size, m
ysize=100000; % Vertical size, m
Nx=35; % Horizontal resolution
Ny=45; % Vertical resolution
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step,m
% Define coordinates of different staggered points
% Basic points
x=0:dx:xsize; % Horizontal coordinates of grid points, m
y=0:dy:ysize; % Vertical coordinates of grid points, m
% Vx-Points
xvx=0:dx:xsize+dx; % horizontal
yvx=-dy/2:dy:ysize+dy/2; % vertical
% Vy-Points
xvy=-dx/2:dx:xsize+dx/2; % horizontal
yvy=0:dy:ysize+dy; % Vertical
% P-Points
xpr=-dx/2:dx:xsize+dx/2; % Horizontal
ypr=-dy/2:dy:ysize+dy/2; % Verical

gy=10; % Vertical gravity, m/s^2
ETA=1e+19; % Viscosity, Pa*s


% Define dinsity in vy-points
RHOvy=zeros(Ny+1,Nx+1); % Density, kg/m^3
for j=1:1:Nx+1
    for i=1:1:Ny+1
        % Define density depending on xvy(j) and yvy(i)
        RHOvy(i,j)=3300; % Medium
        d=((xvy(j)-xsize/2)^2+(yvy(i)-ysize/2)^2)^0.5; % distance to the centre from the given nodal point
        if(d<20000)
            RHOvy(i,j)=3200;
        end
    end
end

% Define global matrixes
N=(Nx+1)*(Ny+1)*3; % Global number of unknown
L=sparse(N,N); % Left parts
R=zeros(N,1); % Right parts

% Composing global matrixes
for j=1:1:Nx+1
    for i=1:1:Ny+1
    % Global indexes
    kvx=((j-1)*(Ny+1)+(i-1))*3+1;
    kvy=kvx+1;
    kpr=kvx+2;
    
    % Equations for Vx
    if(i==1 || i==Ny+1 || j==1 || j==Nx || j==Nx+1)
        % BC: vx=0
        L(kvx,kvx)=1; % Left part
        R(kvx,1)=0; % Right part
        
    else
        % x-Stokes
        % Left part
        L(kvx,kvx-(Ny+1)*3)=ETA/dx^2;  
        L(kvx,kvx-3)=ETA/dy^2;  
        L(kvx,kvx)=-2*ETA/dx^2-2*ETA/dy^2;  
        L(kvx,kvx+3)=ETA/dy^2;  
        L(kvx,kvx+(Ny+1)*3)=ETA/dx^2;  
        L(kvx,kpr)=1/dx;
        L(kvx,kpr+(Ny+1)*3)=-1/dx;
        % Right part
        R(kvx,1)=0;
    end

    % Equations for Vy
    if(i==1 || i==Ny || i==Ny+1 || j==1 || j==Nx+1)
        % BC: vy=0
        L(kvy,kvy)=1; % Left part
        R(kvy,1)=0; % Right part
    else
        % y-Stokes
        % Left part
        L(kvy,kvy-(Ny+1)*3)=ETA/dx^2;  
        L(kvy,kvy-3)=ETA/dy^2;  
        L(kvy,kvy)=-2*ETA/dx^2-2*ETA/dy^2;  
        L(kvy,kvy+3)=ETA/dy^2;  
        L(kvy,kvy+(Ny+1)*3)=ETA/dx^2;  
        L(kvy,kpr)=1/dy;
        L(kvy,kpr+3)=-1/dy;
        % Right part
        R(kvy,1)=-RHOvy(i,j)*gy;
    end

    % Equations for P
    if(i==1 || i==Ny+1 || j==1 || j==Nx+1 || (i==2 && j==2))
        % BC: P=0
        L(kpr,kpr)=1; % Left part
        R(kpr,1)=0; % Right part
        % BC for P in a single point: P=dy/2*gy*3300
        if(i==2 && j==2)
            L(kpr,kpr)=1; % Left part
            R(kpr,1)=0;
        end
        
    else
        % continuity
        % Left part
        L(kpr,kvx-(Ny+1)*3)=-1/dx;  
        L(kpr,kvx)=1/dx;  
        L(kpr,kvy-3)=-1/dy;  
        L(kpr,kvy)=1/dy;  
      % Right part
        R(kpr,1)=0;
    end

    end
end

% Solving global matrixes
S=L\R;

% Reload solutions to vx(), vy(), pr()
vx=zeros(Ny+1,Nx+1);
vy=zeros(Ny+1,Nx+1);
pr=zeros(Ny+1,Nx+1);
for j=1:1:Nx+1
    for i=1:1:Ny+1
    % Global indexes
    kvx=((j-1)*(Ny+1)+(i-1))*3+1;
    kvy=kvx+1;
    kpr=kvx+2;
    % Reload solutions
    vx(i,j)=S(kvx);
    vy(i,j)=S(kvy);
    pr(i,j)=S(kpr);
    end
end

% Visualize RHOv(), vx(), vy(), pr()
% Visualization of geometrical arrays OMEGA(), PSI(), vx(), vy()
figure(1);clf;colormap('Jet')
% Drawing RHO()
subplot(2,2,1)
pcolor(xvy,yvy,RHOvy)
shading flat
colorbar
axis ij image;
title('density, kg/m^3')
% Drawing OMEGA()
subplot(2,2,2)
pcolor(xpr,ypr,pr)
shading interp
colorbar
axis ij image;
title('Pressure, Pa')
% Drawing vx()
subplot(2,2,3)
pcolor(xvx,yvx,vx)
shading interp
colorbar
axis ij image;
title('vx-velocity, m/s')
% Drawing vy()
subplot(2,2,4)
pcolor(xvy,yvy,vy)
shading interp
colorbar
axis ij image;
title('vy-velocity, m/s')
print ('-djpeg', '-r150','staggered');


% ETA=1e+19
aaa(1,1)=pr(7,5);    % 3.74968014815624e+08
aaa(2,1)=vx(7,5);    %-5.95393267966115e-10
aaa(3,1)=vy(7,5);    % 6.24355570975085e-10 










            
            
            
            
            
            
            







