%% Temperature eq in Lagrangian model
% implicit
% Yuan Xie, 23.11.2022
% Homework 9 

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
xT=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of grid points, m
yT=-dy/2:dy:ysize+dy/2; % Vertical coordinates of grid points, m
K=3;
T0=zeros(Ny1,Nx1);
Tdt=zeros(Ny1,Nx1);
RHOCp=zeros(Ny1,Nx1);


% Define density and temperature field
C=[xsize/2,ysize/2];% center point
radius=20*1000; % radium of density anormaly,m

for i=1:1:Ny1
    for j=1:1:Nx1
        x1=xT(j);
        y1=yT(i); % location of each point
        d=sqrt((x1-C(1))^2+(y1-C(2))^2); % distance between center point and each point
        if d <= radius
            RHO=3200;
            Cp=1100;
            RHOCp(i,j)=RHO*Cp;  % density
            T0(i,j)=1773;         
        else
            RHO=3300;
            Cp=1000;            
            RHOCp(i,j)=RHO*Cp;  % density
            T0(i,j)=1573;
        end
    end
end

% 3) Define global matrixes L(), R()
N=Nx1*Ny1; % Global number of unknowns
LT=sparse(N,N); % Matrix of coefficients (left part)
RT=zeros(N,1); % Vector of right parts


Nmax=10; % Number of time steps
for niter=1:1:Nmax
   
    dt=min(dx,dy)^2/(4*K/min(min(RHOCp)));
    % 4) Composing global matrixes L(), R()
    % Going through all points of the grid (2D meant two loops)
    % composing respective equations
    for j=1:1:Nx1
        for i=1:1:Ny1
            g=(j-1)*Ny1+i;
            % BC
            if(i==1)
                LT(g,g)=1;
                LT(g,g+1)=1;
                RT(g,1)=1573*2;
            elseif(i==Ny1)
                LT(g,g)=1;
                LT(g,g-1)=1;
                RT(g,1)=1573*2;
            elseif(j==1)
                LT(g,g)=1;
                LT(g,g+Ny1)=-1;
                RT(g,1)=0;
            elseif(j==Nx1)
                LT(g,g)=1;
                LT(g,g-Ny1)=-1;
                RT(g,1)=0;
            else
                %                    T2
                %                   i-1,j
                %                    |
                %                    |
                %         T1---------T3---------T5
                %         i,j-1     i,j         i,j+1
                %                    |
                %                    |
                %                    T4
                %                   i+1,j
                LT(g,g-Ny1)=-K/dx^2;  %T1
                LT(g,g-1)=-K/dy^2;    %T2
                LT(g,g)=RHOCp(i,j)/dt-K*(-2/dx^2-2/dy^2); %T3
                LT(g,g+1)=-K/dy^2;    %T4
                LT(g,g+Ny1)=-K/dx^2;  %T5
                RT(g,1)=RHOCp(i,j)/dt*T0(i,j);
            end
        end
    end
    
    % 5) Solving matrixes
    S=LT\RT;
    
    % 6) Reload S--> PSI
    for j=1:1:Nx1
        % Second loop - vertical index i
        for i=1:1:Ny1        
            g=(j-1)*Ny1+i;
            % Reload solution
            Tdt(i,j)=S(g);
        end
    end
    T0=Tdt;    
end




figure(1); clf;
% Colormap
colormap('Jet');
subplot(1,2,1)
pcolor(xT,yT,RHOCp)
colorbar
shading flat;
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down

subplot(1,2,2)
pcolor(xT,yT,T0)
colorbar
shading interp
axis ij image; % Image sizes propotional to coordinates, vertical axis upside down


aaa(1,1)=T0(17,15);
format long e
disp(aaa)
% aaa(1,1)=TDT(17,15); %[1701.50831679871], K=3, 2D 100x100km, 35x45 explicit
% aaa(1,1)=TDT(17,15); %[1704.92596509904], K=3, 2D 100x100km, 35x45  implicit


