%% Solving 1D advection equation
% Mark in cell
% using finite differences in time
% Homework 6
% Yuan Xie, 03/11/2022

% 1) Clear memory and figures
clear all
clf

% 2) Define model
xsize=100; % Horizontal model size, m

% 2.1) Define Eulerian grid
Nx=101; % Model resolution
dx=xsize/(Nx-1); % Grid size, m
x=0:dx:xsize; % Coordinates of Eularian (immobile) grid points, m

% 2.2) Define Lagrangian grid markers
Nxm=1001; % Number of Lagrangian grid points
dxm=xsize/(Nxm-1); % Number of Lagrangian grid points
xm0=0:dxm:xsize; % Initial coordinates of Lagrangian grid points

% 2.3) Define initial density in Lagrangian points
RHOm=zeros(1,Nxm);
for m=1:1:Nxm
    RHOm(m)=3200; % background
    if(xm0(m)>0.4*xsize && xm0(m)<0.6*xsize)
        RHOm(m)=3300;
    end
end

% 2.4)  define timestep
dt=0.5;
% Coordinates of Lagrangian points for the next moment of time
xmdt=xm0;
ntimesteps=20; % Number of time steps

% Time step loop
for t=1:1:ntimesteps
    
    % Interplolation RHO from makers to Eularian nodes(weight average) 
    RHO0=zeros(1,Nx);
    for j=1:1:Nx
        k=0;
        Wtmj=0;
        rho=0;
        for m=1:1:Nxm
            dXmj=abs(x(j)-xm0(m));
            if dXmj <= dx
                k=k+1;
                Wtmj(k)=1-dXmj/dx;
                rho(k)=RHOm(m)*Wtmj(k);
            end
        end
        RHO0(j)=sum(rho)/sum(Wtmj);
    end

    % Define advection velocity
    % velocity in nodes
    vx=zeros(1,Nx);
    for j=1:1:Nx
        vx(j)=1+sin((x(j)-x(1))/xsize*pi*5)*0.1*(RHO0(j)-3250)/50;
    end
    
    % Interpolate velocity from nodes to makers
    vxm=zeros(1,Nxm);
    vxm(end)=vx(end);
    for m=1:1:Nxm-1
        j=fix(xm0(m)/dx)+1;
        dXmj=x(j)-xm0(m);
        Wtmj=1-abs(dXmj)/dx;
        Wtmj1=1-Wtmj;
        vxm(m)=Wtmj*vx(j)+Wtmj1*vx(j+1);
    end

    % Going through all lagrangian points
    for m=1:1:Nxm
        xmdt(m)=xm0(m)+vxm(m)*dt;
        % Periodic BC
        if(xmdt(m)>xsize)
            xmdt(m)=xmdt(m)-xsize;
        end
        if(xmdt(m)<0)
            xmdt(m)=xmdt(m)+xsize;
        end
    end
    
    
    %     Plot results
    figure(1)
    subplot(2,1,1)
    plot(xm0,RHOm,'o r')
    hold on
    plot(x,RHO0,'-o b')
    hold off
    axis([0 xsize 3150 3350])
    
    subplot(2,1,2)
    plot(x,vx,'-o b')
    axis([0 xsize 0.9 1.1])
    pause(1)
    
    % going to the next time steps
    xm0=xmdt;
end


aaa(1,1)=RHO0(51); % 3285.09563976021
aaa(2,1)=vx(51);  % 1.07019127952042
format long e
disp(aaa(1))
disp('correct answer: 3285.09563976021')
disp(aaa(2))
disp('correct answer: 1.07019127952042')
