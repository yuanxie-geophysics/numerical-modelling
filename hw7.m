%% Solving stocks and continuity equation
% using FD with staggered grid
% markers in cell 2D

% 1) Clear memory and figures
clear all;
clf

% 2) Define numerical model
% 2.1) define in staggerd grid
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


% 2.2) define lagrangian grid makers
Nxm=(Nx-1)*5; % Horizontal grid resolution
Nym=(Ny-1)*5; % Vertical grid resolution
Nm=Nxm*Nym;  % Number of markers
RHOm=zeros(Nm,1); % density in makers
ETAm=zeros(Nm,1); % viscosity in makers
xm=zeros(Nm,1); % Horizontal coordinates
ym=zeros(Nm,1); % Vertical coordinates
dxm=xsize/Nxm; % horizontal step, m
dym=ysize/Nym; % vertical step, m

% 2.3) define density on markers
C=[xsize/2,ysize/2];% center point
radius=20*1000; % radius of density anormaly,m
m=1;
for jm=1:1:Nxm
    for im=1:1:Nym
        % marker coordinates
        %         xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        %         ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        xm(m)=dxm/2+(jm-1)*dxm;
        ym(m)=dym/2+(im-1)*dym;
        % distance
        d=sqrt((xm(m)-C(1))^2+(ym(m)-C(2))^2); % distance between center point and each point
        if d <= radius
            RHOm(m)=3200;  % density
            ETAm(m)=1e18;
        else
            RHOm(m)=3300;  % density
            ETAm(m)=1e19;
        end
        m=m+1;
    end
end

% 2.4) Define global matrixes L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts


% 2.5) Boundary conditions: free slip=-1; No Slip=1
bc=-1;

% 2.6) Define time step
ntimesteps=10; % Number of time steps

% 3) Time step loop
for t=1:1:ntimesteps
    
    % 3.1) interpolate rho, eta, etaP from makers
    rho=zeros(Ny1,Nx1); % Density, kg/m^3, on Vy nodes
    eta=zeros(Ny1,Nx1); % Viscosity, Pa*s, on basic nodes
    etaP=zeros(Ny1,Nx1); % Viscosity, Pa*s, on pressure nodes
    
    RHOwtsum =zeros(Ny,Nx);
    WtRHOsum =zeros(Ny,Nx);
    ETAwtsum =zeros(Ny,Nx);
    WtETAsum =zeros(Ny,Nx);
    ETAPwtsum=zeros(Ny,Nx);
    WtETAPsum=zeros(Ny,Nx);
    % calculate the weight
    for m=1:1:Nm
        % for rho
        j=fix((xm(m)-xVy(1))/dx)+1;
        i=fix((ym(m)-yVy(1))/dy)+1;
        if i>(Ny-1)
            i=Ny-1;
        end
        if j>(Nx-1)
            j=Nx-1;
        end
        % distance
        dxmij=xm(m)-xVy(j);
        dymij=ym(m)-yVy(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        %  i,j
        RHOwtsum(i,j)=RHOwtsum(i,j)+RHOm(m)*wtmij;
        WtRHOsum(i,j)=WtRHOsum(i,j)+wtmij;
        %  i,j+1
        RHOwtsum(i,j+1)=RHOwtsum(i,j+1)+RHOm(m)*wtmij1;
        WtRHOsum(i,j+1)=WtRHOsum(i,j+1)+wtmij1;
        %  i+1,j
        RHOwtsum(i+1,j)=RHOwtsum(i+1,j)+RHOm(m)*wtmi1j;
        WtRHOsum(i+1,j)=WtRHOsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        RHOwtsum(i+1,j+1)=RHOwtsum(i+1,j+1)+RHOm(m)*wtmi1j1;
        WtRHOsum(i+1,j+1)=WtRHOsum(i+1,j+1)+wtmi1j1;
        
        % for eta
        j=fix((xm(m)-x(1))/dx)+1;
        i=fix((ym(m)-y(1))/dy)+1;
        if i>(Ny-1)
            i=Ny-1;
        end
        if j>(Nx-1)
            j=Nx-1;
        end
        % distance
        dxmij=xm(m)-x(j);
        dymij=ym(m)-y(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        %  i,j
        ETAwtsum(i,j)=ETAwtsum(i,j)+ETAm(m)*wtmij;
        WtETAsum(i,j)=WtETAsum(i,j)+wtmij;
        %  i,j+1
        ETAwtsum(i,j+1)=ETAwtsum(i,j+1)+ETAm(m)*wtmij1;
        WtETAsum(i,j+1)=WtETAsum(i,j+1)+wtmij1;
        %  i+1,j
        ETAwtsum(i+1,j)=ETAwtsum(i+1,j)+ETAm(m)*wtmi1j;
        WtETAsum(i+1,j)=WtETAsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        ETAwtsum(i+1,j+1)=ETAwtsum(i+1,j+1)+ETAm(m)*wtmi1j1;
        WtETAsum(i+1,j+1)=WtETAsum(i+1,j+1)+wtmi1j1;
        
        % for etaP
        j=fix((xm(m)-xP(1))/dx)+1;
        i=fix((ym(m)-yP(1))/dy)+1;
        if i>(Ny-1)
            i=Ny-1;
        end
        if j>(Nx-1)
            j=Nx-1;
        end
        % distance
        dxmij=xm(m)-xP(j);
        dymij=ym(m)-yP(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        %  i,j
        ETAPwtsum(i,j)=ETAPwtsum(i,j)+ETAm(m)*wtmij;
        WtETAPsum(i,j)=WtETAPsum(i,j)+wtmij;
        %  i,j+1
        ETAPwtsum(i,j+1)=ETAPwtsum(i,j+1)+ETAm(m)*wtmij1;
        WtETAPsum(i,j+1)=WtETAPsum(i,j+1)+wtmij1;
        %  i+1,j
        ETAPwtsum(i+1,j)=ETAPwtsum(i+1,j)+ETAm(m)*wtmi1j;
        WtETAPsum(i+1,j)=WtETAPsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        ETAPwtsum(i+1,j+1)=ETAPwtsum(i+1,j+1)+ETAm(m)*wtmi1j1;
        WtETAPsum(i+1,j+1)=WtETAPsum(i+1,j+1)+wtmi1j1;
    end
    % Calculate rho, eta, etaP
    for j=1:Nx
        for i=1:Ny
            rho(i,j)=RHOwtsum(i,j)/WtRHOsum(i,j);
            eta(i,j)=ETAwtsum(i,j)/WtETAsum(i,j);
            etaP(i,j)=ETAPwtsum(i,j)/WtETAPsum(i,j);
        end
    end

    % 3.2) Going through all points of the grid
    % composing respective equations
    for j=1:1:Nx1
        for i=1:1:Ny1
            % Define global index gVx, gVy, gP
            gVx=((j-1)*Ny1+(i-1))*3+1; % Vx
            gVy=gVx+1; % Vy
            gP =gVx+2; % P
            
            % 3.2a) Eq. for Vx (ghost BC x-stokes)
            if j==Nx1 % ghost nodes
                L(gVx,gVx)=1; % Left part
                R(gVx,1)=0; % Right part
            else
                if (i==1 || j==1 || j==Nx || i==Ny1 )    % Boundary Condition
                    L(gVx,gVx)=1; % Left part
                    R(gVx,1)=0; % Right part
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
            % 3.2b) Eq. for Vy (ghost BC x-stokes)
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
            
            % 3.2c) Eq. for P (ghost BC Continuity-eq)
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
    % 3.3) Solving matrixes
    S=L\R;
    % 3.4) Reload S--> vx,vy,P
    Vx=zeros(Ny1,Nx1); % Vx, m/s
    Vy=zeros(Ny1,Nx1); % Vy, m/s
    P=zeros(Ny1,Nx1);
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
    % 3.5) plot
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
    pause(0.01)
    
    % 3.6) define timestep
    Vxmax=max(max(abs(Vx)));
    Vymax=max(max(abs(Vy)));
    dt=1e+30;
    if (dt>0.5*dx/Vxmax)
        dt=0.5*dx/Vxmax;
    end
    if (dt>0.5*dy/Vymax)
        dt=0.5*dy/Vymax;
    end
    
    % 3.7) move markers
    % 3.7a) interpolate velocity to makers
    for m=1:1:Nm
        % for vx points
        j=fix((xm(m)-xVx(1))/dx)+1;
        i=fix((ym(m)-yVx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % distance
        dxmij=xm(m)-xVx(j);
        dymij=ym(m)-yVx(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        % vx velocity
        vxm=Vx(i,j)*wtmij+Vx(i+1,j)*wtmi1j+...
            Vx(i,j+1)*wtmij1+Vx(i+1,j+1)*wtmi1j1;
        
        % for vy points
        j=fix((xm(m)-xVy(1))/dx)+1;
        i=fix((ym(m)-yVy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        dxmij=xm(m)-xVy(j);
        dymij=ym(m)-yVy(i);
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wti1j1=(dxmij/dx)*(dymij/dy);
        % vy velocity
        vym=Vy(i,j)*wtmij+Vy(i+1,j)*wtmi1j+...
            Vy(i,j+1)*wtmij1+Vy(i+1,j+1)*wtmi1j1;
        
        % move makers
        xm(m)=xm(m)+dt*vxm;
        ym(m)=ym(m)+dt*vym;
    end   
end

% FREE SLIP 
aaa(1,1)=P(7,5);   % STEP1:  1374925092.91932        STEP10: 1374890082.55988
aaa(2,1)=Vx(7,5);  % STEP1:  -1.73770057946995e-09   STEP10: -2.06597098103458e-09
aaa(3,1)=Vy(7,5);  % STEP1: 1.78377334721631e-09     STEP10: 2.20165981984881e-09

format long e
disp(aaa)
