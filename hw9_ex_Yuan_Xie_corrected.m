%% Temperature eq in Lagrangian model
% explicit
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

Nmax=10; % Number of time steps
for niter=1:1:Nmax
    for j=1:1:Nx1
        for i=1:1:Ny1
            dt=min(dx,dy)^2/(4*K/min(min(RHOCp)));
%             % ERROR !!!: BC should be done after all internal points
%             % and use Tdt not T0
%             % BC
%             if(i==1)
%                  Tdt(i,j)=1573*2-T0(i+1,j);
%             elseif(i==Ny1)
%                 Tdt(i,j)=1573*2-T0(i-1,j);
%             elseif(j==1)
%                 Tdt(i,j)=T0(i,j+1);
%             elseif ( j==Nx1 )
%                 Tdt(i,j)=T0(i,j-1);
%             else
             if(j>1 && j<Nx+1 && i>1 && i<Ny+1)
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
                T1=T0(i,j-1);
                T2=T0(i-1,j);
                T3=T0(i,j);
                T4=T0(i+1,j);
                T5=T0(i,j+1);
                Tdt(i,j)=K*dt/RHOCp(i,j)*((T5-2*T3+T1)/dx^2+(T4-2*T3+T2)/dy^2)+T3;
            end
        end
    end
    % ERROR !!!: BC should be done after all internal points
    % and use Tdt not T0
    for j=1:1:Nx1
        for i=1:1:Ny1
            dt=min(dx,dy)^2/(4*K/min(min(RHOCp)));
            % BC
            if(i==1)
                 Tdt(i,j)=1573*2-Tdt(i+1,j);
            elseif(i==Ny1)
                Tdt(i,j)=1573*2-Tdt(i-1,j);
            elseif(j==1)
                Tdt(i,j)=Tdt(i,j+1);
            elseif ( j==Nx1 )
                Tdt(i,j)=Tdt(i,j-1);
            end
        end
    end

    
    T0=Tdt;
end

% reload 
T0=Tdt;


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
aaa(1,1)=Tdt(17,15); %[1701.50831679871], K=3, 2D 100x100km, 35x45 explicit
% aaa(1,1)=TDT(17,15); %[1704.92596509904], K=3, 2D 100x100km, 35x45  implicit


