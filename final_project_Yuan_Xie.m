%% Numerical Modelling --Final Project
% Yuan Xie
% 22/01/2023


% 1) Clear memory and figures
clear all;
clf

% 2) Define numerical model
% constant values
gy = 10; % Gravity acceleration, m/s^2

% 2.1) define in staggerd grid
xsize = 500*1000; % Horizontal model size, m
ysize = 400*1000; % Vertical model size, m
Nx = 101; % Horizontal grid resolution
Ny = 81; % Vertical grid resolution
Nx1 = Nx+1;
Ny1 = Ny+1;

% coordinates
dx = xsize/(Nx-1); % Horizontal grid step, m
dy = ysize/(Ny-1); % Vertical grid step, m
x = 0:dx:xsize+dx; % Horizontal coordinates of basic grid points, m
y = 0:dy:ysize+dy; % Vertical coordinates of basic grid points, m
xVx = 0:dx:xsize+dx; % Horizontal coordinates of vx grid points, m
yVx = -dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
xVy = -dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yVy = 0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
xP = -dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yP = -dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m

% for Stokes equations
Vx = zeros(Ny1,Nx1); % Vx, m/s
Vy = zeros(Ny1,Nx1); % Vy, m/s
P = zeros(Ny1,Nx1); % Pressure, Pa
rhoVy = zeros(Ny1,Nx1); % Density, kg/m^3, on Vy nodes
eta = zeros(Ny1,Nx1); % Viscosity, Pa*s, on basic nodes
etaP = zeros(Ny1,Nx1); % Viscosity, Pa*s, on pressure nodes
% for thermal eq
rhoT = zeros(Ny1,Nx1); % Density, kg/m^3, on T nodes
T0 = zeros(Ny1,Nx1);
Tdt = zeros(Ny1,Nx1);
HR = zeros(Ny1,Nx1);
Hs = zeros(Ny1,Nx1);
Ha = zeros(Ny1,Nx1);
Alpha = zeros(Ny1,Nx1);
kx = zeros(Ny1,Nx1);
ky = zeros(Ny1,Nx1);
kp = zeros(Ny1,Nx1);
RHOCP = zeros(Ny1,Nx1);
% for temperature subgrid diffusion
NRD = zeros(Ny1, Nx1);
DTsubgrid = zeros(Ny1, Nx1);
DTremaining = zeros(Ny1, Nx1);
% for stress and strain rate
Exx = zeros(Ny1, Nx1);
Sxx = zeros(Ny1, Nx1);
Eyy = zeros(Ny1, Nx1);
Syy = zeros(Ny1, Nx1);
Exy = zeros(Ny1, Nx1);
Sxy = zeros(Ny1, Nx1);

% 2.2) define lagrangian grid makers
Nxm = (Nx-1)*5; % Horizontal grid resolution
Nym = (Ny-1)*5; % Vertical grid resolution
Nm = Nxm*Nym;  % Number of markers
dxm = xsize/Nxm; % horizontal step, m
dym = ysize/Nym; % vertical step, m
TYm = zeros(Nm,1);
% for continuity eq
Pm = zeros(Nm,1);
RHOm = zeros(Nm,1); % density in makers
ETAm = zeros(Nm,1); % viscosity in makers
xm = zeros(Nm,1); % Horizontal coordinates
ym = zeros(Nm,1); % Vertical coordinates
%for thermal eq
Tm = zeros(Nm,1);
RHOCPm = zeros(Nm,1);
Km = zeros(Nm,1);
ALPHAm = zeros(Nm,1); % Thermal expansion
HRm = zeros(Nm,1);  % Radioative
% for temperature subgrid diffusion
Tmn = zeros(Nm, 1);
dTm0 = zeros(Nm, 1);
dTmdt = zeros(Nm, 1);
NRDm = zeros(Nm, 1);
DTmsubgrid = zeros(Nm, 1);
DTmremaining = zeros(Nm, 1);
% for stress and strain rate
Exxm = zeros(Nm,1);
Eyym = zeros(Nm,1);
Exym = zeros(Nm,1);
EIIm = zeros(Nm,1);

% 2.3) define properties on markers
% material properties
    % rock type 
    % [1] air
    % [2] dry olivine for Necking area
    % [3] dry olivine for slab and plate
    % [4] wet olivine for mantle
%          1            2           3         4 
rho0    = [  1        3400       3400     3300     ];
Hr      = [  0        2e-8       2e-8     3e-8     ];
strength= [  0        2e+7       1e+8     5e+7     ];
eta0    = [  1e+18    1e+23      1e+23    1e+20    ];
cp      = [  3.3e+6   1000       1000     1000     ];
alpha   = [  0        3e-5       3e-5     3e-5     ];
betta   = [  0        1e-11      1e-11    1e-11    ];
AD      = [  0        2.5e-17    2.5e-17  2e-21    ];
n_rock  = [  0        3.5        3.5      4        ];
Va      = [  0        8e-6       8e-6     4e-6     ];
Ea      = [  0        532000     532000   471000   ];
% Limits for ETA
ETAmax = 1e24;
ETAmin = 1e18;
% initial model setup
m=1;
for jm=1:1:Nxm
    for im=1:1:Nym
        % marker coordinates
        %         xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        %         ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        xm(m)=dxm/2+(jm-1)*dxm;
        ym(m)=dym/2+(im-1)*dym;
        % Asthenosphere  rock type [4]
        Tm(m)=1573;
        %         RHOm(m)=rho0(4);
        RHOm(m)=3250;
        ETAm(m)=eta0(4);  % Viscosity
        Cp=cp(4);
        RHOCPm(m)=Cp*RHOm(m);
        Km(m)=0.73+1293/(Tm(m)+77);
        ALPHAm(m)=alpha(4);
        HRm(m)=Hr(4);
        TYm(m)=4;
        % Sticky air layer rock type [1]
        if (ym(m) < 50*1000 )
            RHOm(m)=rho0(1);
            ETAm(m)=eta0(1);
            Cp=cp(1);
            RHOCPm(m)=Cp*RHOm(m);
            Tm(m)=273;
            Km(m)=3000;
            ALPHAm(m)=alpha(1);
            HRm(m)=Hr(1);
            TYm(m)=1;
        end
        
        % lithosphere layer  Plate rock type [3]
        depth_top=50; % km
        depth_bottom=100; %km
        if ( ym(m) >= depth_top*1000 && ym(m)<= depth_bottom*1000)
            RHOm(m)=rho0(3);
            ETAm(m)=eta0(3);
            Cp=cp(3);
            RHOCPm(m)=Cp*RHOm(m);
            ALPHAm(m)=alpha(3);
            Tm(m)=(1573-273)/((depth_bottom-depth_top)*1000)*...
                (ym(m)-depth_top*1000)+273;
            Km(m)=0.73+1293/(Tm(m)+77);
            HRm(m)=Hr(3);
            TYm(m)=3;
            % 
        end
        % slab  rock type [3]
        slab_bottom=50+50+150;
            if ( ym(m) >= depth_top*1000 && ym(m)<=100*1000 && ...
                    xm(m) <= (230*1000+ ym(m)-50*1000) && ...
                    xm(m) >= (220*1000+(ym(m)-50*1000)/5))
                RHOm(m)=rho0(3);
                ETAm(m)=eta0(3);
                Cp=cp(3);
                RHOCPm(m)=Cp*RHOm(m);
                ALPHAm(m)=alpha(3);
                HRm(m)=Hr(3);
                Tm(m)=(1573-273)/60/1000*(60*1000-(xm(m)-220*1000-(ym(m)-100*1000)))+273;
                Km(m)=0.73+1293/(Tm(m)+77);
                TYm(m)=3;
            end
            % for necking area rock type [2]
            if (  ym(m)>100*1000 && ym(m)<= 150*1000 && ...
                    xm(m) >= (220*1000+(ym(m)-100*1000)) && ...
                    xm(m) <= (230*1000+ ym(m)-50*1000))
                RHOm(m)=rho0(2);
                ETAm(m)=eta0(2);
                Cp=cp(2);
                RHOCPm(m)=Cp*RHOm(m);
                ALPHAm(m)=alpha(2);
                HRm(m)=Hr(2);
                Tm(m)=(1573-273)/60/1000*(60*1000-(xm(m)-220*1000-(ym(m)-100*1000)))+273;
                Km(m)=0.73+1293/(Tm(m)+77);
                TYm(m)=2;
            end
            % for slab in asthenosphere rock type [3]
            if  (xm(m) >= (220*1000+(ym(m)-100*1000)) && ...
                    xm(m) <= (230*1000+ ym(m)-50*1000) && ...
                    ym(m)>150*1000 && ym(m)<= slab_bottom*1000 )
                RHOm(m)=rho0(3);
                ETAm(m)=eta0(3);
                Cp=cp(3);
                RHOCPm(m)=Cp*RHOm(m);
                ALPHAm(m)=alpha(3);
                HRm(m)=Hr(3);
                Tm(m)=(1573-273)/60/1000*(60*1000-(xm(m)-220*1000-(ym(m)-100*1000)))+273;
                Km(m)=0.73+1293/(Tm(m)+77);
                TYm(m)=3;
            end
            m=m+1;
    end
end

% 3) Define global matrixes L(), R()
N = Nx1*Ny1*3; % Global number of unknowns
L = sparse(N,N); % Matrix of coefficients (left part)
R = zeros(N,1); % Vector of right parts
NT = Nx1*Ny1; % Global number of unknowns
LT = sparse(NT,NT); % Matrix of coefficients (left part)
RT = zeros(NT,1); % Vector of right parts

% 3.1) Boundary conditions: free slip=-1; No Slip=1
bc = -1;
dispmax = 0.5;
Tchangemax = 50; % Max temp`erature change per time step

% 3.2) Define time step
ntimesteps=50; % Number of time steps
dt=1e+11;
timesum=0;

% 4) Time step loop
for t=1:1:ntimesteps
    % 4.1) interpolate rho, eta, etaP from makers to nodes
    %  Clear transport properties arrays
    % for stokes eq.
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
    KPwtsum=zeros(Ny1,Nx1);
    WtKPsum=zeros(Ny1,Nx1);
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
        j=fix((xm(m)-xVx(1))/dx)+1;
        i=fix((ym(m)-yVx(1))/dy)+1;
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
        dxmij=xm(m)-xVx(j);
        dymij=ym(m)-yVx(i);
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
        j=fix((xm(m)-xVy(1))/dx)+1;
        i=fix((ym(m)-yVy(1))/dy)+1;
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
        dxmij=xm(m)-xVy(j);
        dymij=ym(m)-yVy(i);
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
        
        % for RHOCp
        %  i,j
        RHOCPwtsum(i,j)=RHOCPwtsum(i,j)+RHOCPm(m)*wtmij;
        WtRHOCPsum(i,j)=WtRHOCPsum(i,j)+wtmij;
        %  i,j+1
        RHOCPwtsum(i,j+1)=RHOCPwtsum(i,j+1)+RHOCPm(m)*wtmij1;
        WtRHOCPsum(i,j+1)=WtRHOCPsum(i,j+1)+wtmij1;
        %  i+1,j
        RHOCPwtsum(i+1,j)=RHOCPwtsum(i+1,j)+RHOCPm(m)*wtmi1j;
        WtRHOCPsum(i+1,j)=WtRHOCPsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        RHOCPwtsum(i+1,j+1)=RHOCPwtsum(i+1,j+1)+RHOCPm(m)*wtmi1j1;
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
        TmRhocpwtsum(i,j)=TmRhocpwtsum(i,j)+Tm(m)*RHOCPm(m)*wtmij;
        Rhocpwtsum(i,j)=Rhocpwtsum(i,j)+RHOCPm(m)*wtmij;
        %  i,j+1
        TmRhocpwtsum(i,j+1)=TmRhocpwtsum(i,j+1)+Tm(m)*RHOCPm(m)*wtmij1;
        Rhocpwtsum(i,j+1)=Rhocpwtsum(i,j+1)+RHOCPm(m)*wtmij1;
        %  i+1,j
        TmRhocpwtsum(i+1,j)=TmRhocpwtsum(i+1,j)+Tm(m)*RHOCPm(m)*wtmi1j;
        Rhocpwtsum(i+1,j)=Rhocpwtsum(i+1,j)+RHOCPm(m)*wtmi1j;
        %  i+1,j+1
        TmRhocpwtsum(i+1,j+1)=TmRhocpwtsum(i+1,j+1)+Tm(m)*RHOCPm(m)*wtmi1j1;
        Rhocpwtsum(i+1,j+1)=Rhocpwtsum(i+1,j+1)+RHOCPm(m)*wtmi1j1;
        
        % for Hr
        %  i,j
        Hrwtsum(i,j)=Hrwtsum(i,j)+HRm(m)*wtmij;
        WtHrsum(i,j)=WtHrsum(i,j)+wtmij;
        %  i,j+1
        Hrwtsum(i,j+1)=Hrwtsum(i,j+1)+HRm(m)*wtmij1;
        WtHrsum(i,j+1)=WtHrsum(i,j+1)+wtmij1;
        %  i+1,j
        Hrwtsum(i+1,j)=Hrwtsum(i+1,j)+HRm(m)*wtmi1j;
        WtHrsum(i+1,j)=WtHrsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        Hrwtsum(i+1,j+1)=Hrwtsum(i+1,j+1)+HRm(m)*wtmi1j1;
        WtHrsum(i+1,j+1)=WtHrsum(i+1,j+1)+wtmi1j1;
        
        % for Alpha
        %  i,j
        ALPHAwtsum(i,j)=ALPHAwtsum(i,j)+ALPHAm(m)*wtmij;
        WtALPHAsum(i,j)=WtALPHAsum(i,j)+wtmij;
        %  i,j+1
        ALPHAwtsum(i,j+1)=ALPHAwtsum(i,j+1)+ALPHAm(m)*wtmij1;
        WtALPHAsum(i,j+1)=WtALPHAsum(i,j+1)+wtmij1;
        %  i+1,j
        ALPHAwtsum(i+1,j)=ALPHAwtsum(i+1,j)+ALPHAm(m)*wtmi1j;
        WtALPHAsum(i+1,j)=WtALPHAsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        ALPHAwtsum(i+1,j+1)=ALPHAwtsum(i+1,j+1)+ALPHAm(m)*wtmi1j1;
        WtALPHAsum(i+1,j+1)=WtALPHAsum(i+1,j+1)+wtmi1j1;
        
        % for Kp
        KPwtsum(i,j)=KPwtsum(i,j)+Km(m)*wtmij;
        WtKPsum(i,j)=WtKPsum(i,j)+wtmij;
        %  i,j+1
        KPwtsum(i,j+1)=KPwtsum(i,j+1)+Km(m)*wtmij1;
        WtKPsum(i,j+1)=WtKPsum(i,j+1)+wtmij1;
        %  i+1,j
        KPwtsum(i+1,j)=KPwtsum(i+1,j)+Km(m)*wtmi1j;
        WtKPsum(i+1,j)=WtKPsum(i+1,j)+wtmi1j;
        %  i+1,j+1
        KPwtsum(i+1,j+1)=KPwtsum(i+1,j+1)+Km(m)*wtmi1j1;
        WtKPsum(i+1,j+1)=WtKPsum(i+1,j+1)+wtmi1j1;
    end
    % Calculate rho, eta, etaP
    for j=1:Nx1
        for i=1:Ny1
            if (WtRHOsum(i,j)>0)
                rhoVy(i,j)=RHOwtsum(i,j)/WtRHOsum(i,j);
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
                RHOCP(i,j)=RHOCPwtsum(i,j)/WtRHOCPsum(i,j);
            end
            if (WtRHOTsum(i,j)>0)
                rhoT(i,j)=RHOTwtsum(i,j)/WtRHOTsum(i,j);
            end
            if (Rhocpwtsum(i,j)>0)
                T0(i,j)=TmRhocpwtsum(i,j)/Rhocpwtsum(i,j);
            end
            if (WtHrsum(i,j)>0)
                HR(i,j)=Hrwtsum(i,j)/WtHrsum(i,j);
            end
            if (WtALPHAsum(i,j)>0)
                Alpha(i,j)=ALPHAwtsum(i,j)/WtALPHAsum(i,j);
            end
            if (WtKPsum(i,j)>0)
                kp(i,j)=KPwtsum(i,j)/WtKPsum(i,j);
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
        % 4.2) Going through all points of the grid
        % composing respective equations
        for j=1:1:Nx1
            for i=1:1:Ny1
                % Define global index gVx, gVy, gP
                gVx=((j-1)*Ny1+(i-1))*3+1; % Vx
                gVy=gVx+1; % Vy
                gP =gVx+2; % P
                % 4.2a) Eq. for Vx (ghost BC x-stokes)
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
                % 4.2b) Eq. for Vy (ghost BC x-stokes)
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
                        rho1=rhoVy(i,j-1);
                        rho2=rhoVy(i-1,j);
                        rho4=rhoVy(i+1,j);
                        rho5=rhoVy(i,j+1);
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
                        R(gVy,1)=-rhoVy(i,j)*gy;
                    end
                end
                
                % 4.2c) Eq. for P (ghost BC Continuity-eq)
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
        
        
        % 4.3) Solving matrixes
        S=L\R;
        
        % 4.4) Reload S--> vx,vy,P
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
        
        
        % 4.5) define timestep (not in final iteration)
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
                Hs(i,j)=Sxx(i,j)*Exx(i,j)+Syy(i,j)*Eyy(i,j)+...
                    2*(Exy(i,j)*Sxy(i,j)+Exy(i-1,j)*Sxy(i-1,j)+...
                    Exy(i,j-1)*Sxy(i,j-1)+Exy(i-1,j-1)*Sxy(i-1,j-1))/4;             
            end
        end

        
        % Calculate Ha 
        % on T nodes
        for j=1:1:Nx1
            for i=2:1:Ny1
                vy=(Vy(i,j)+Vy(i-1,j))/2;
                Ha(i,j)=Alpha(i,j)*rhoT(i,j)*gy*vy;
            end
        end   
        
        % 5) thermal solver
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
                    LT(gT,gT)=RHOCP(i,j)/dt+((kx2+kx1)/dx^2+...
                        (ky2+ky1)/dy^2); %T3
                    LT(gT,gT+1)=-ky2/dy^2;    %T4
                    LT(gT,gT+Ny1)=-kx2/dx^2;  %T5
                    % right side
                    RT(gT,1)=RHOCP(i,j)/dt*T0(i,j)+HR(i,j)+Hs(i,j)+T0(i,j)/2*Ha(i,j);
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
    
   
    % 6) plot
    figure(1);colormap('Jet');
    subplot(3,4,1)
    pcolor(xVy,yVy,rhoVy)
    hold on
    [c,h]=contour(xVy,yVy(1:Ny),rhoVy(1:Ny,:),[-1 1650.5 5000],'m','LineWidth',2);
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
    pcolor(xP,yP,T0)
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
    pcolor(xP,yP,RHOCP)
    colorbar
    shading flat;
    axis ij image; % Image sizes propotional to coordinates, vertical axis upside down
    title('RhoCp')
    
    subplot(3,4,8)
    pcolor(xVx,yVx,log10(kx))
    colorbar
    shading flat;
    axis ij image;
    title('kx')
    
    subplot(3,4,9)
    pcolor(xVy,yVy,Alpha)
    colorbar
    shading flat;
    axis ij image;
    title('ALPHA')
    
    subplot(3,4,10)
    pcolor(xP,yP,HR)
    colorbar
    shading flat;
    axis ij image;
    title('HR')
    
    subplot(3,4,11)
    pcolor(xP,yP,Hs)
    colorbar
    shading flat;
    axis ij image;
    title('Hs')
    
    subplot(3,4,12)
    pcolor(xP,yP,Ha.*(Tdt+T0)/2)
    colorbar
    shading flat;
    axis ij image;
    title('Ha')
    
    
    % 3.7a) interpolate temperature to makers
    if (t == 1 )
        for m=1:1:Nm
            % for T points
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
            
            Tm(m)=Tdt(i,j)*wtmij+Tdt(i+1,j)*wtmi1j+...
                Tdt(i,j+1)*wtmij1+Tdt(i+1,j+1)*wtmi1j1;
        end
    end
    
    
    DTMSUBGRIDPSUM = zeros(Ny1, Nx1);
    if(t > 1)
        % a. Calculate NRD on P/T nodes
        for j=1:1:Nx1
            for i=1:1:Ny1
                NRD(i,j) = exp((-kp(i,j)/RHOCP(i,j))*((1/dx^2)+(1/dy^2)*dt));
            end
        end   
         % b. DTmsubgrid on makers    
        for m=1:1:Nm
            j=fix((xm(m)-xP(1))/dx)+1;
            i=fix((ym(m)-yP(1))/dy)+1;
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
            dxmij=xm(m)-xP(j);
            dymij=ym(m)-yP(i);
            % weight
            wtmij = (1-(dxmij/dx))*(1-(dymij/dy));
            wtmi1j = (1-(dxmij/dx))*(dymij/dy);
            wtmij1 = (dxmij/dx)*(1-(dymij/dy));
            wtmi1j1 = (dxmij/dx)*(dymij/dy);
            
            % Interpolation
            % c. Uppdate Tm(m) with DTm & d. Interpolate DTm to nodes
            Tmn(m) = T0(i,j)*wtmij + T0(i+1,j)*wtmi1j + T0(i,j+1)*wtmij1 + ...
                T0(i+1,j+1)*wtmi1j1;
            dTm0(m) = Tm(m)-Tmn(m);
            NRDm(m) =  NRD(i,j)*wtmij + NRD(i+1,j)*wtmi1j + ...
                NRD(i,j+1)*wtmij1 + NRD(i+1,j+1)*wtmi1j1;
            dTmdt(m) = dTm0(m)* NRDm(m);
            dTmdt(m) = dTm0(m)* NRDm(m);
            DTmsubgrid(m) = dTmdt(m)-dTm0(m);
            
            % Interpolate DTmsubgrid to P nodes
            DTMSUBGRIDPSUM(i,j) = DTMSUBGRIDPSUM(i,j) + RHOCPm(m)*DTmsubgrid(m)*wtmij;
            DTMSUBGRIDPSUM(i+1,j) = DTMSUBGRIDPSUM(i+1,j) + RHOCPm(m)*DTmsubgrid(m)*wtmi1j;
            DTMSUBGRIDPSUM(i,j+1) = DTMSUBGRIDPSUM(i,j+1) + RHOCPm(m)*DTmsubgrid(m)*wtmij1;
            DTMSUBGRIDPSUM(i+1,j+1) = DTMSUBGRIDPSUM(i+1, j+1) + RHOCPm(m)*DTmsubgrid(m)*wtmi1j1;
        end
        
        % Calculate DTsubgrid(i,j) in P/T nodes
        for j=1:1:Nx+1
            for i = 1:1:Ny+1
                if(Rhocpwtsum(i,j)>0)
                    DTsubgrid(i,j) = DTMSUBGRIDPSUM(i,j)/Rhocpwtsum(i,j);
                end
            end
        end
        
        % Calculate DTremaining
        DTremaining = dT - DTsubgrid;
        
        % Interpolate DTmremaining
        for m=1:1:Nm
            if(xm(m)>0 && xm(m)<xsize && ym(m)>0 && ym(m)<ysize)
                % Interpolation of Tmn
                j = fix((xm(m)-xP(1))/dx)+1;
                i = fix((ym(m)-yP(1))/dy)+1;
                
                % Distance to upper left velocity point
                dxmij = xm(m)-xP(j);
                dymij = ym(m)-yP(i);
                
                % Distance based weight
                wtmij = (1-(dxmij/dx))*(1-(dymij/dy));
                wtmi1j = (1-(dxmij/dx))*(dymij/dy);
                wtmij1 = (dxmij/dx)*(1-(dymij/dy));
                wtmi1j1 = (dxmij/dx)*(dymij/dy);
                
                %Interpolation of DTmremaining from pressure nodule points
                DTmremaining(m) = DTremaining(i,j)*wtmij + ...
                    DTremaining(i+1,j)*wtmi1j + DTremaining(i,j+1)*wtmij1 +...
                    DTremaining(i+1,j+1)*wtmi1j1;
                
                Tm(m) = Tm(m) + DTmsubgrid(m) + DTmremaining(m);
            end
        end
    end

    
    % Interpolate Exx Eyy Exy P to makers
    for m=1:1:Nm
        j=fix((xm(m)-xP(1))/dx)+1;
        i=fix((ym(m)-yP(1))/dy)+1;
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
        dxmij=xm(m)-xP(j);
        dymij=ym(m)-yP(i);
        
        % Distance based weight
        wtmij = (1-(dxmij/dx))*(1-(dymij/dy));
        wtmi1j = (1-(dxmij/dx))*(dymij/dy);
        wtmij1 = (dxmij/dx)*(1-(dymij/dy));
        wtmi1j1 = (dxmij/dx)*(dymij/dy);
        
        %Interpolation
        Pm(m) = P(i,j)*wtmij + P(i+1,j)*wtmi1j + P(i,j+1)*wtmij1 + P(i+1,j+1)*wtmi1j1;
        % Exxm and Eyym from P nodes
        Exxm(m) = Exx(i,j)*wtmij + Exx(i+1,j)*wtmi1j + Exx(i,j+1)*wtmij1 + Exx(i+1,j+1)*wtmi1j1;
        Eyym(m) = Eyy(i,j)*wtmij + Eyy(i+1,j)*wtmi1j + Eyy(i,j+1)*wtmij1 + Eyy(i+1,j+1)*wtmi1j1;
        
        % Interpolate pressure, Exy from basic nodes
        j = fix((xm(m)-x(1))/dx)+1;
        i = fix((ym(m)-y(1))/dy)+1;
        
        dxmij = xm(m)-x(j);
        dymij = ym(m)-y(i);
        
        % weight
        wtmij = (1-(dxmij/dx))*(1-(dymij/dy));
        wtmi1j = (1-(dxmij/dx))*(dymij/dy);
        wtmij1 = (dxmij/dx)*(1-(dymij/dy));
        wtmi1j1 = (dxmij/dx)*(dymij/dy);
        
        %Interpolation
        Exym(m) = Exy(i,j)*wtmij + Exy(i+1,j)*wtmi1j + Exy(i,j+1)*wtmij1 + Exy(i+1,j+1)*wtmi1j1;
        
        % Calculate second strain rate invariant
        EIIm(m) = sqrt(0.5*(Exxm(m)^2 + Eyym(m)^2) + Exym(m)^2);
        
        % update properties
        if(TYm(m) > 1) 
            RHOm(m) = rho0(TYm(m))*(1+betta(TYm(m))*(Pm(m)-1e+5))/(1+alpha(TYm(m))*(Tm(m)-273));
            RHOCPm(m) = RHOm(m)*cp(TYm(m));
            Km(m) = 0.73 + (1293/(Tm(m) + 77));
            ETAm(m) = 0.5 * (1/((AD(TYm(m)))^(1/n_rock(TYm(m))))) * (EIIm(m)^((1/n_rock(TYm(m)))-1)) * exp((Ea(TYm(m))+Pm(m)*Va(TYm(m)))/(n_rock(TYm(m))*8.314*Tm(m)));
            if (ETAm(m) > strength(TYm(m))/(2*EIIm(m)) )
                ETAm(m) = strength(TYm(m))/(2*EIIm(m));
            end
            if(ETAm(m) > ETAmax)
                ETAm(m) = ETAmax;
            elseif(ETAm(m) < ETAmin)
                ETAm(m) = ETAmin;
            end
        end     
    end
      
    % 3.7) move markers
    % 3.7a) interpolate velocity to makers
    for m=1:1:Nm
               
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
    t
    timesum=timesum+dt;
end