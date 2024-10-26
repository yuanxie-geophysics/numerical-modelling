%% Solving stocks and continuity equation
% using FD with staggered grid
% markers in cell 2D
% Sticky air Free surface
% Runge Kutta
% with thermal mechanical
% Yuan Xie, 15/12/2022

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

%for thermal eq
xT=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of grid points, m
yT=-dy/2:dy:ysize+dy/2; % Vertical coordinates of grid points, m
T0=zeros(Ny1,Nx1);
Tdt=zeros(Ny1,Nx1);
Hr=zeros(Ny1,Nx1);
Hs=zeros(Ny1,Nx1);
Ha=zeros(Ny1,Nx1);
alpha=zeros(Ny1,Nx1);
xKx=xVx;
yKx=yVx;
xKy=xVy;
yKy=yVy;

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

%for thermal eq
Tm=zeros(Nm,1);
RHOCpm=zeros(Nm,1);
Km=zeros(Nm,1);
alpham=zeros(Nm,1); % Thermal expansion
Hrm=zeros(Nm,1);  % Radioative

% 2.3) define properties on markers
C=[xsize/2,ysize*0.7];% center point of the plume
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
        if d <= radius   % INSIDE the plume
            RHOm(m)=3200;  % density
            ETAm(m)=1e18;  % Viscosity
            Cp=1100;
            RHOCpm(m)=Cp*RHOm(m);
            Km(m)=2;
            Tm(m)=1773;
            alpham(m)=2.2*10^-5;
            Hrm(m)=3*10^-8;
        else  % in Asthenosphere
            RHOm(m)=3300;  % density
            ETAm(m)=1e19;  % Viscosity
            Cp=1000;
            RHOCpm(m)=Cp*RHOm(m);
            Km(m)=3;
            Tm(m)=1573;
            alpham(m)=2.1*10^-5;
            Hrm(m)=2*10^-8;
        end
        % Sticky air layer
        if (ym(m)<ysize*0.2)
            RHOm(m)=1;
            ETAm(m)=1e+17;
            RHOCpm(m)=3.3e+6;
%             RHOCpm(m)=3e+6;
            Km(m)=3000;
            Tm(m)=273;
            alpham(m)=0;
            Hrm(m)=0;
        end
        % lithosphere layer
        if (ym(m)>=20*1000 && ym(m)<=40*1000)
            RHOm(m)=3350;
            ETAm(m)=1e21;
            Cp=900;
            RHOCpm(m)=Cp*RHOm(m);
            Km(m)=4;
            alpham(m)=2*10^-5;
            % calcalate temperature in m
            Tm(m)=(1573-273)/(20*1000)*(ym(m)-20*1000)+273;
            Hrm(m)=10^-8;
        end
        m=m+1;
    end
end

% 2.4) Define global matrixes L(), R()
N=Nx1*Ny1*3; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts
NT=Nx1*Ny1; % Global number of unknowns
LT=sparse(NT,NT); % Matrix of coefficients (left part)
RT=zeros(NT,1); % Vector of right parts

% 2.5) Boundary conditions: free slip=-1; No Slip=1
bc=-1;
dispmax=0.5;
Tchangemax=50; % Max temperature change per time step
% 2.6) Define time step
ntimesteps=10; % Number of time steps
dt=1e+11;
timesum=0;
% 3) Time step loop
for t=1:1:ntimesteps
    % 3.1) interpolate rho, eta, etaP from makers to nodes
    %  Clear transport properties arrays
    % for stokes eq.
    rho=zeros(Ny1,Nx1); % Density, kg/m^3, on Vy nodes
    rhoT=zeros(Ny1,Nx1); % Density, kg/m^3, on T nodes
    eta=zeros(Ny1,Nx1); % Viscosity, Pa*s, on basic nodes
    etaP=zeros(Ny1,Nx1); % Viscosity, Pa*s, on pressure nodes
    
    % for thermal eq.
    Dt=zeros(Ny1,Nx1);
    kx=zeros(Ny1,Nx1);
    ky=zeros(Ny1,Nx1);
    RHOCp=zeros(Ny1,Nx1);
    Hr=zeros(Ny1,Nx1);
    % Clear wights for basic nodes
    RHOwtsum =zeros(Ny1,Nx1);
    WtRHOsum =zeros(Ny1,Nx1);
    ETAwtsum =zeros(Ny1,Nx1);
    WtETAsum =zeros(Ny1,Nx1);
    ETAPwtsum=zeros(Ny1,Nx1);
    WtETAPsum=zeros(Ny1,Nx1);
    KXwtsum=zeros(Ny1,Nx1);
    WtKXsum=zeros(Ny1,Nx1);
    KYwtsum=zeros(Ny1,Nx1);
    WtKYsum=zeros(Ny1,Nx1);
    RHOCPwtsum=zeros(Ny1,Nx1);
    WtRHOCPsum=zeros(Ny1,Nx1);
    RHOTwtsum=zeros(Ny1,Nx1);
    WtRHOTsum=zeros(Ny1,Nx1);
    TmRhocpwtsum=zeros(Ny1,Nx1);
    Rhocpwtsum=zeros(Ny1,Nx1);
    Hrwtsum=zeros(Ny1,Nx1);
    WtHrsum=zeros(Ny1,Nx1);
    ALPHAwtsum=zeros(Ny1,Nx1);
    WtALPHAsum=zeros(Ny1,Nx1);
    
    % calculate the weight
    for m=1:1:Nm
        % for rho
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
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
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
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
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
        
        % for Kx
        j=fix((xm(m)-xKx(1))/dx)+1;
        i=fix((ym(m)-yKx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % distance
        dxmij=xm(m)-xKx(j);
        dymij=ym(m)-yKx(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        %  i,j
        KXwtsum(i,j)=KXwtsum(i,j)+Km(m)*wtmij;
        WtKXsum(i,j)=WtKXsum(i,j)+wtmij;
        %  i,j+1
        KXwtsum(i,j+1)=KXwtsum(i,j+1)+Km(m)*wtmij1;
        WtKXsum(i,j+1)=WtKXsum(i,j+1)+wtmij1;
        %  i+1,j
        KXwtsum(i+1,j)=KXwtsum(i+1,j)+Km(m)*wtmi1j;
        WtKXsum(i+1,j)=WtKXsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        KXwtsum(i+1,j+1)=KXwtsum(i+1,j+1)+Km(m)*wtmi1j1;
        WtKXsum(i+1,j+1)=WtKXsum(i+1,j+1)+wtmi1j1;
        
        % for KY
        j=fix((xm(m)-xKy(1))/dx)+1;
        i=fix((ym(m)-yKy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % distance
        dxmij=xm(m)-xKy(j);
        dymij=ym(m)-yKy(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        %  i,j
        KYwtsum(i,j)=KYwtsum(i,j)+Km(m)*wtmij;
        WtKYsum(i,j)=WtKYsum(i,j)+wtmij;
        %  i,j+1
        KYwtsum(i,j+1)=KYwtsum(i,j+1)+Km(m)*wtmij1;
        WtKYsum(i,j+1)=WtKYsum(i,j+1)+wtmij1;
        %  i+1,j
        KYwtsum(i+1,j)=KYwtsum(i+1,j)+Km(m)*wtmi1j;
        WtKYsum(i+1,j)=WtKYsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        KYwtsum(i+1,j+1)=KYwtsum(i+1,j+1)+Km(m)*wtmi1j1;
        WtKYsum(i+1,j+1)=WtKYsum(i+1,j+1)+wtmi1j1;
        
        % for T0 & RHOCp & Hr & rhoT
        j=fix((xm(m)-xT(1))/dx)+1;
        i=fix((ym(m)-yT(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % distance
        dxmij=xm(m)-xT(j);
        dymij=ym(m)-yT(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        
        % for RHOCp
        %  i,j
        RHOCPwtsum(i,j)=RHOCPwtsum(i,j)+RHOCpm(m)*wtmij;
        WtRHOCPsum(i,j)=WtRHOCPsum(i,j)+wtmij;
        %  i,j+1
        RHOCPwtsum(i,j+1)=RHOCPwtsum(i,j+1)+RHOCpm(m)*wtmij1;
        WtRHOCPsum(i,j+1)=WtRHOCPsum(i,j+1)+wtmij1;
        %  i+1,j
        RHOCPwtsum(i+1,j)=RHOCPwtsum(i+1,j)+RHOCpm(m)*wtmi1j;
        WtRHOCPsum(i+1,j)=WtRHOCPsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        RHOCPwtsum(i+1,j+1)=RHOCPwtsum(i+1,j+1)+RHOCpm(m)*wtmi1j1;
        WtRHOCPsum(i+1,j+1)=WtRHOCPsum(i+1,j+1)+wtmi1j1;
        
        % for rhoT
        %  i,j
        RHOTwtsum(i,j)=RHOTwtsum(i,j)+RHOm(m)*wtmij;
        WtRHOTsum(i,j)=WtRHOTsum(i,j)+wtmij;
        %  i,j+1
        RHOTwtsum(i,j+1)=RHOTwtsum(i,j+1)+RHOm(m)*wtmij1;
        WtRHOTsum(i,j+1)=WtRHOTsum(i,j+1)+wtmij1;
        %  i+1,j
        RHOTwtsum(i+1,j)=RHOTwtsum(i+1,j)+RHOm(m)*wtmi1j;
        WtRHOTsum(i+1,j)=WtRHOTsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        RHOTwtsum(i+1,j+1)=RHOTwtsum(i+1,j+1)+RHOm(m)*wtmi1j1;
        WtRHOTsum(i+1,j+1)=WtRHOTsum(i+1,j+1)+wtmi1j1;
        
        % for T0
        %  i,j
        TmRhocpwtsum(i,j)=TmRhocpwtsum(i,j)+Tm(m)*RHOCpm(m)*wtmij;
        Rhocpwtsum(i,j)=Rhocpwtsum(i,j)+RHOCpm(m)*wtmij;
        %  i,j+1
        TmRhocpwtsum(i,j+1)=TmRhocpwtsum(i,j+1)+Tm(m)*RHOCpm(m)*wtmij1;
        Rhocpwtsum(i,j+1)=Rhocpwtsum(i,j+1)+RHOCpm(m)*wtmij1;
        %  i+1,j
        TmRhocpwtsum(i+1,j)=TmRhocpwtsum(i+1,j)+Tm(m)*RHOCpm(m)*wtmi1j;
        Rhocpwtsum(i+1,j)=Rhocpwtsum(i+1,j)+RHOCpm(m)*wtmi1j;
        %  i+1,j+1
        TmRhocpwtsum(i+1,j+1)=TmRhocpwtsum(i+1,j+1)+Tm(m)*RHOCpm(m)*wtmi1j1;
        Rhocpwtsum(i+1,j+1)=Rhocpwtsum(i+1,j+1)+RHOCpm(m)*wtmi1j1;
        
        % for Hr
        %  i,j
        Hrwtsum(i,j)=Hrwtsum(i,j)+Hrm(m)*wtmij;
        WtHrsum(i,j)=WtHrsum(i,j)+wtmij;
        %  i,j+1
        Hrwtsum(i,j+1)=Hrwtsum(i,j+1)+Hrm(m)*wtmij1;
        WtHrsum(i,j+1)=WtHrsum(i,j+1)+wtmij1;
        %  i+1,j
        Hrwtsum(i+1,j)=Hrwtsum(i+1,j)+Hrm(m)*wtmi1j;
        WtHrsum(i+1,j)=WtHrsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        Hrwtsum(i+1,j+1)=Hrwtsum(i+1,j+1)+Hrm(m)*wtmi1j1;
        WtHrsum(i+1,j+1)=WtHrsum(i+1,j+1)+wtmi1j1;
        
        % for Alpha
        %  i,j
        ALPHAwtsum(i,j)=ALPHAwtsum(i,j)+alpham(m)*wtmij;
        WtALPHAsum(i,j)=WtALPHAsum(i,j)+wtmij;
        %  i,j+1
        ALPHAwtsum(i,j+1)=ALPHAwtsum(i,j+1)+alpham(m)*wtmij1;
        WtALPHAsum(i,j+1)=WtALPHAsum(i,j+1)+wtmij1;
        %  i+1,j
        ALPHAwtsum(i+1,j)=ALPHAwtsum(i+1,j)+alpham(m)*wtmi1j;
        WtALPHAsum(i+1,j)=WtALPHAsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        ALPHAwtsum(i+1,j+1)=ALPHAwtsum(i+1,j+1)+alpham(m)*wtmi1j1;
        WtALPHAsum(i+1,j+1)=WtALPHAsum(i+1,j+1)+wtmi1j1;
    end
    % Calculate rho, eta, etaP
    for j=1:Nx1
        for i=1:Ny1
            if (WtRHOsum(i,j)>0)
                rho(i,j)=RHOwtsum(i,j)/WtRHOsum(i,j);
            end
            if (WtETAsum(i,j)>0)
                eta(i,j)=ETAwtsum(i,j)/WtETAsum(i,j);
            end
            if (WtETAPsum(i,j)>0)
                etaP(i,j)=ETAPwtsum(i,j)/WtETAPsum(i,j);
            end
            if (WtKXsum(i,j)>0)
                kx(i,j)=KXwtsum(i,j)/WtKXsum(i,j);
            end
            if (WtKYsum(i,j)>0)
                ky(i,j)=KYwtsum(i,j)/WtKYsum(i,j);
            end
            if (WtRHOCPsum(i,j)>0)
                RHOCp(i,j)=RHOCPwtsum(i,j)/WtRHOCPsum(i,j);
            end
            if (WtRHOTsum(i,j)>0)
                rhoT(i,j)=RHOTwtsum(i,j)/WtRHOTsum(i,j);
            end
            if (Rhocpwtsum(i,j)>0)
                T0(i,j)=TmRhocpwtsum(i,j)/Rhocpwtsum(i,j);
            end
            if (WtHrsum(i,j)>0)
                Hr(i,j)=Hrwtsum(i,j)/WtHrsum(i,j);
            end
            if (WtALPHAsum(i,j)>0)
                alpha(i,j)=ALPHAwtsum(i,j)/WtALPHAsum(i,j);
            end
        end
    end
    % BC fot T0
    T0(2:Ny,1)=T0(2:Ny,2); % Left BC: dT/dx=0
    T0(2:Ny,Nx1)=T0(2:Ny,Nx); % Right BC: dT/dx=0
    T0(1,:)=2*273-T0(2,:); % Top BC: T=273
    T0(Ny1,:)=2*1573-T0(Ny,:); % Bottom BC: T=1573
    
    
    % Introducing thermomechanical iteration
    nitermax=2;
    % % probing increase of the time step
    dt=dt*1.1; % coefficient 1.1 implies gradual (by 10% max) increase in the time step size
    for niter=1:1:nitermax
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
                        % density
                        rho1=rho(i,j-1);
                        rho2=rho(i-1,j);
                        rho4=rho(i+1,j);
                        rho5=rho(i,j+1);
                        drhox=(rho5-rho1)/2/dx;
                        drhoy=(rho4-rho2)/2/dy;
                        % Left part
                        L(gVy,gVy-Ny1*3)=eta1/dx^2; % vy1
                        L(gVy,gVy-3)=2*etaP1/dy^2; % vy2
                        L(gVy,gVy)=-2*(etaP2+etaP1)/dy^2-...
                            (eta2+eta1)/dx^2-...
                            dt*drhoy*gy; % vy3
                        L(gVy,gVy+3)=2*etaP2/dy^2;  % vy4
                        L(gVy,gVy+Ny1*3)=eta2/dx^2;  % vy5
                        L(gVy,gVx-Ny1*3)=eta1/dx/dy-dt*drhox*gy/4; % vx1
                        L(gVy,gVx-Ny1*3+3)=-eta1/dx/dy-dt*drhox*gy/4; % vx2
                        L(gVy,gVx)=-eta2/dx/dy-dt*drhox*gy/4;  % vx3
                        L(gVy,gVx+3)=eta2/dx/dy-dt*drhox*gy/4; % vx4
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
                        R(gP,1)=1e5; % Right part
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
        
        
        % 3.6) define timestep (not in final iteration)
        if (niter < nitermax)
            Vxmax=max(max(abs(Vx)));
            Vymax=max(max(abs(Vy)));
%             dt = 1e+30;
            if (dt>dispmax*dx/Vxmax)
                dt=dispmax*dx/Vxmax;
                dtmechanicalvx=dt
            end
            if (dt>0.5*dy/Vymax)
                dt=0.5*dy/Vymax;
                dtmechanicalvy=dt
            end
        end
        
        % calculate Hs
        %  
        Sxx=zeros(Ny1,Nx1); % define on T nodes
        Syy=zeros(Ny1,Nx1); % define on T nodes
        Exx=zeros(Ny1,Nx1); % define on T nodes
        Eyy=zeros(Ny1,Nx1); % define on T nodes
        
        for j=2:1:Nx1
            for i=2:1:Ny1
                % exx=dVx/dx  
                Exx(i,j)=(Vx(i,j)-Vx(i,j-1))/dx;
                % Sxx=2*etaP*exx
                Sxx(i,j)=2*etaP(i,j)*Exx(i,j);
                % eyy=dVy/dy
                Eyy(i,j)=(Vy(i,j)-Vy(i-1,j))/dy;
                % Syy=2*etaP*eyy
                Syy(i,j)=2*etaP(i,j)*Eyy(i,j);                 
            end
        end

        Sxy=zeros(Ny1,Nx1); % define on basic nodes
        Exy=zeros(Ny1,Nx1); % define on basic nodes        
        for j=1:1:Nx
            for i=1:1:Ny
                % exy=1/2*(dVx/dy+dVy/dx)
                Exy(i,j)=1/2*((Vx(i+1,j)-Vx(i,j))/dy+...
                    (Vy(i,j+1)-Vy(i,j))/dx);
                % Sxy=2*eta*exy;
                Sxy(i,j)=2*eta(i,j)*Exy(i,j);
            end
        end
        
        for j=2:1:Nx1
            for i=2:1:Ny1
                % ERROR !!!: average product needed
%                 EXY=(Exy(i,j)+Exy(i-1,j)+Exy(i,j-1)+Exy(i-1,j-1))/4;
%                 SXY=(Sxy(i,j)+Sxy(i-1,j)+Sxy(i,j-1)+Sxy(i-1,j-1))/4;
%                 Hs(i,j)=Sxx(i,j)*Exx(i,j)+Syy(i,j)*Eyy(i,j)+...
%                     2*SXY*EXY;             
                Hs(i,j)=Sxx(i,j)*Exx(i,j)+Syy(i,j)*Eyy(i,j)+...
                    2*(Exy(i,j)*Sxy(i,j)+Exy(i-1,j)*Sxy(i-1,j)+Exy(i,j-1)*Sxy(i,j-1)+Exy(i-1,j-1)*Sxy(i-1,j-1))/4;             
            end
        end

        
        % Calculate Ha 
        % on T nodes
        for j=1:1:Nx1
            for i=2:1:Ny1
                vy=(Vy(i,j)+Vy(i-1,j))/2;
                % ERROR !!!: should be implicit approach
%                 if t==1
%                     Ha(i,j)=T0(i,j)*alpha(i,j)*rhoT(i,j)*gy*vy;
%                 else
%                     Ha(i,j)=(T0(i,j)+Tdt(i,j))/2*alpha(i,j)*rhoT(i,j)*gy*vy;
%                 end
                 Ha(i,j)=alpha(i,j)*rhoT(i,j)*gy*vy;
            end
        end   
        
        % 3.7) thermal solver
        for j=1:1:Nx1
            for i=1:1:Ny1
                gT=(j-1)*Ny1+i; % T                
                % 3.2d) Eq. for thermal equations
                % BC
                if(i==1)
                    LT(gT,gT)=1;
                    LT(gT,gT+1)=1;
                    RT(gT,1)=273*2;
                elseif(i==Ny1)
                    LT(gT,gT)=1;
                    LT(gT,gT-1)=1;
                    RT(gT,1)=1573*2;
                elseif(j==1)
                    LT(gT,gT)=1;
                    LT(gT,gT+Ny1)=-1;
                    RT(gT,1)=0;
                elseif(j==Nx1)
                    LT(gT,gT)=1;
                    LT(gT,gT-Ny1)=-1;
                    RT(gT,1)=0;
                else
                    % left side
                    kx1=kx(i,j-1);
                    kx2=kx(i,j);
                    ky1=ky(i-1,j);
                    ky2=ky(i,j);
                    LT(gT,gT-Ny1)=-kx1/dx^2;  %T1
                    LT(gT,gT-1)=-ky1/dy^2;    %T2
                    LT(gT,gT)=RHOCp(i,j)/dt+((kx2+kx1)/dx^2+...
                        (ky2+ky1)/dy^2); %T3
                    LT(gT,gT+1)=-ky2/dy^2;    %T4
                    LT(gT,gT+Ny1)=-kx2/dx^2;  %T5
                    % right side
                % ERROR !!!: should be implicit approach
%                    RT(gT,1)=RHOCp(i,j)/dt*T0(i,j)+Hr(i,j)+Hs(i,j)+Ha(i,j);
                   RT(gT,1)=RHOCp(i,j)/dt*T0(i,j)+Hr(i,j)+Hs(i,j)+T0(i,j)/2*Ha(i,j);
                   LT(gT,gT)=LT(gT,gT)-1/2*Ha(i,j);
                end
            end
        end
        ST=LT\RT;
        Tdt=zeros(Ny1,Nx1);
        for j=1:1:Nx1
            % Second loop - vertical index i
            for i=1:1:Ny1
                % Define global index g
                gT=(j-1)*Ny1+i;
                Tdt(i,j)=ST(gT);
            end
        end
        dT=Tdt-T0;
        
        % Do not change time step during the last iteration
        if (niter<nitermax)
            DTmax=max(max(abs(dT)));
            if (DTmax>Tchangemax)
                dt=dt*Tchangemax/DTmax/1.1; % for safety we make dt 10% less than the allowed dt
                dtthermal=dt;
            end
        end
    end
    niter
    dt
    % Exit thermomechanical iteration
    % 3.5) plot
    figure(1);colormap('Jet');
    subplot(3,4,1)
    pcolor(xVy,yVy,rho)
    hold on
    [c,h]=contour(xVy,yVy(1:Ny),rho(1:Ny,:),[-1 1650.5 5000],'m','LineWidth',2);
    clabel(c,h,'LabelSpacing',1100,'FontSize',10,'Color','m');
    hold off
    shading flat
    colorbar
    axis ij image;
    title('Density, kg/m^3')
    
    
    subplot(3,4,2)
    pcolor(x,y,log10(eta))
    shading flat;
    axis ij image;
    colorbar
    title('log ETAB(Pa.s)')
    
    subplot(3,4,3)
    pcolor(xT,yT,T0)
    shading interp;
    axis ij image;
    colorbar
    title('Temperature(K)')
    
    subplot(3,4,4)
    pcolor(xP,yP,P)
    shading interp;
    axis ij image;
    colorbar
    title('Pressure(Pa)')
    
    subplot(3,4,5)
    pcolor(xVx,yVx,Vx)
    shading interp;
    axis ij image;
    colorbar
    title('vx(m/s)')
    
    subplot(3,4,6)
    pcolor(xVy,yVy,Vy)
    shading interp;
    axis ij image;
    colorbar
    title('vy(m/s)')
    pause(0.01)
    
    subplot(3,4,7)
    pcolor(xT,yT,RHOCp)
    colorbar
    shading flat;
    axis ij image; % Image sizes propotional to coordinates, vertical axis upside down
    title('RhoCp')
    
    subplot(3,4,8)
    pcolor(xKx,yKx,log10(kx))
    colorbar
    shading flat;
    axis ij image;
    title('kx')
    
    subplot(3,4,9)
    pcolor(xKy,yKy,alpha)
    colorbar
    shading flat;
    axis ij image;
    title('ALPHA')
    
    subplot(3,4,10)
    pcolor(xT,yT,Hr)
    colorbar
    shading flat;
    axis ij image;
    title('Hr')
    
    subplot(3,4,11)
    % ERROR !!!: wrong coordinates
%     pcolor(xKy,yKy,Hs)
    pcolor(xT,yT,Hs)
    colorbar
    shading flat;
    axis ij image;
    title('Hs')
    
    subplot(3,4,12)
    % ERROR !!!: implicit approach needed, wrong coordinates
%     pcolor(xKy,yKy,Ha)
    pcolor(xT,yT,Ha.*(Tdt+T0)/2)
    colorbar
    shading flat;
    axis ij image;
    title('Ha')
    
    
    % 3.7) move markers
    % 3.7a) interpolate velocity and temperature to makers
    for m=1:1:Nm
        % for T points
        j=fix((xm(m)-xT(1))/dx)+1;
        i=fix((ym(m)-yT(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % distance
        dxmij=xm(m)-xT(j);
        dymij=ym(m)-yT(i);
        % weights
        wtmij =(1-dxmij/dx)*(1-dymij/dy);
        wtmi1j=(1-dxmij/dx)*(dymij/dy);
        wtmij1=(dxmij/dx)*(1-dymij/dy);
        wtmi1j1=(dxmij/dx)*(dymij/dy);
        % temperature change
        dTm=dT(i,j)*wtmij+dT(i+1,j)*wtmi1j+...
            dT(i,j+1)*wtmij1+dT(i+1,j+1)*wtmi1j1;
        if t==1
            Tm(m)=Tdt(i,j)*wtmij+Tdt(i+1,j)*wtmi1j+...
                Tdt(i,j+1)*wtmij1+Tdt(i+1,j+1)*wtmi1j1;
        else
            Tm(m)=Tm(m)+dTm;
        end
        
        % for vx
        % for RK
        xA=xm(m);
        yA=ym(m);
        % runge kutta
        Vxm=zeros(4,1);
        Vym=zeros(4,1);
        for rk=1:1:4
            % for vx points
            j=fix((xm(m)-xVx(1))/dx)+1;
            i=fix((ym(m)-yVx(1))/dy)+1;
            % distance
            dxmij=xm(m)-xVx(j);
            dymij=ym(m)-yVx(i);
            % weights
            wtmij =(1-dxmij/dx)*(1-dymij/dy);
            wtmi1j=(1-dxmij/dx)*(dymij/dy);
            wtmij1=(dxmij/dx)*(1-dymij/dy);
            wtmi1j1=(dxmij/dx)*(dymij/dy);
            % correction
            correction=0;
            if (dxmij > dx/2 && j<Nx-1)
                correction=1/2*(dxmij/dx-0.5)^2*(((Vx(i,j)-2*Vx(i,j+1)+...
                    Vx(i,j+2))*(1-dymij/dy))+(Vx(i+1,j)-2*Vx(i+1,j+1)+...
                    Vx(i+1,j+2))*dymij/dy);
            elseif (dxmij <= dx/2 && j>1)
                correction=1/2*(dxmij/dx-0.5)^2*(((Vx(i,j-1)-2*Vx(i,j)+...
                    Vx(i,j+1))*(1-dymij/dy))+(Vx(i+1,j-1)-2*Vx(i+1,j)+...
                    Vx(i+1,j+1))*dymij/dy);
            end
            
            % vx velocity
            Vxm(rk)=Vx(i,j)*wtmij+Vx(i+1,j)*wtmi1j+...
                Vx(i,j+1)*wtmij1+Vx(i+1,j+1)*wtmi1j1+correction;
            
            % for vy points
            j=fix((xm(m)-xVy(1))/dx)+1;
            i=fix((ym(m)-yVy(1))/dy)+1;
            dxmij=xm(m)-xVy(j);
            dymij=ym(m)-yVy(i);
            wtmij =(1-dxmij/dx)*(1-dymij/dy);
            wtmi1j=(1-dxmij/dx)*(dymij/dy);
            wtmij1=(dxmij/dx)*(1-dymij/dy);
            wtmi1j1=(dxmij/dx)*(dymij/dy);
            % correction
            correction=0;
            if (dymij > dy/2 && i<Ny-1)
                correction=1/2*(dymij/dy-0.5)^2*(((Vy(i,j)-2*Vy(i+1,j)+...
                    Vy(i+2,j))*(1-dxmij/dx))+(Vy(i,j+1)-2*Vy(i+1,j+1)+...
                    Vy(i+2,j+1))*dxmij/dx);
            elseif (dymij <= dy/2 && i>1)
                correction=1/2*(dymij/dy-0.5)^2*(((Vy(i-1,j)-2*Vy(i,j)+...
                    Vy(i+1,j))*(1-dxmij/dx))+(Vy(i-1,j+1)-2*Vy(i,j+1)+...
                    Vy(i+1,j+1))*dxmij/dx);
            end
            % vy velocity
            Vym(rk)=Vy(i,j)*wtmij+Vy(i+1,j)*wtmi1j+...
                Vy(i,j+1)*wtmij1+Vy(i+1,j+1)*wtmi1j1+correction;
            
            % coordinates of B,C,D points
            if (rk==1 || rk==2)
                xm(m)=xA+dt/2*Vxm(rk);
                ym(m)=yA+dt/2*Vym(rk);
            elseif (rk==3)
                xm(m)=xA+dt*Vxm(rk);
                ym(m)=yA+dt*Vym(rk);
            end
            
        end
        vxmeff=1/6*(Vxm(1)+2*Vxm(2)+2*(Vxm(3))+Vxm(4));
        vymeff=1/6*(Vym(1)+2*Vym(2)+2*(Vym(3))+Vym(4));
        xm(m)=xA;
        ym(m)=yA;
        % move makers
        xm(m)=xm(m)+dt*vxmeff;
        ym(m)=ym(m)+dt*vymeff;
    end
    timesum=timesum+dt;
end

Ha=Ha.*(T0+Tdt)/2;
% FREE SLIP
aaa(1,1)=P(27,12); % 1264697641.62341
aaa(2,1)=Vx(27,12); % -4.58472978899238e-09
aaa(3,1)=Vy(27,12); %  2.53836272700191e-10
aaa(4,1)=Tdt(27,12); % 1574.83958476236
aaa(5,1)=dt;       % 66351281694.0421
aaa(6,1)=timesum; % 66351281694.0421
aaa(7,1)=Hr(27,12);% 2.00480000000000e-08
aaa(8,1)=Hs(27,12);% 4.61750912647993e-07
aaa(9,1)=Ha(27,12);% 2.14943490257245e-07


format long e
disp(aaa)
