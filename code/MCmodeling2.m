%Question 2
close all;
clc
%Kwabena Gyasi Bawuah
%101048814
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%electron spec
 global C

    addpath ../geom2d/geom2d

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                    % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665; %metres (32.1740 ft) per s²
    
    T = 300;
    k = 1.38e-23;
    mn = 0.26*C.m_0; %effective mass
    tmn = 0.2e-12;    % Mean time between collisions
   
    vth = sqrt((2*C.kb*T)/mn);% Thermal velocity
    
    freepath = vth*tmn   % mean free path
    
    ConductorL = 180e-9;
    ConductorW = 80e-9;
    
    dpoints = 5e4;
    ecount = 15; %the number of electron to show on plot
    
    detaT= ConductorW/vth/100;
    sims = 1000;
    
%     Xpos = rand(1,ecount).*ConductorW;
%     Ypos = rand(1,ecount).*ConductorL;
    traj=zeros(sims,ecount*2);
    temp=zeros(sims, 1);
    
    temp(:,1)= 300;
     %------------------------------------
    %Collisions with Mean Free Path
   
    Pscat = 1-exp(-detaT/tmn);
    %to make a probability distribution with mu=0 and sigma=vth
    ProbDistr = makedist('Normal','mu', 0, 'sigma', sqrt(C.kb*T/mn));
    for i = 1: dpoints
        state(i,:)= [ConductorL*rand ConductorW*rand random(ProbDistr) random(ProbDistr)];
    end
   
    %from part 1
    for i = 1 :sims
    state(:,1:2)=state(:,1:2)+detaT.*state(:,3:4);
    %specifying the particles reactions at boundary
    out = state(:,1)> ConductorL;
    state(out,1) = state(out,1)-ConductorL;
    
    out = state(:,2) < 0;
    state(out,2) = -state(out,2);
    state(out,4) = -state(out,4);
    
    out = state(:,2)> ConductorW;
    state(out,2)= 2 * ConductorW - state(out,2);
    state(out,4)= -state(out,4);
    
    out = state(:,1)< 0;
    state(out,1)=state(out,1)+ ConductorL;
    
    out = rand(ecount,1) < Pscat;
    state(out,3:4)=random(ProbDistr,[sum(out),2]);
    
    %varying temp 
    temp(i)=(sum(state(:,3).^2) + sum(state(:,4).^2)).*mn/k/2/dpoints;
    %iterating over array of visible electrons
    for out = 1:ecount 
        traj(i, (2 * out):(2 * out + 1)) = state(out, 1:2);
    end
    %plot on every 5 iteration steps
    if mod(i,5)==0
        figure(3);
        subplot(3,1,1);
        part = sqrt(state(:,3).^2 + state(:,4).^2);
        xlim([0 7e5]);
        ylim([0 2000]);
        histogram(part);
        xlabel('v(m/s)');
        ylabel('Particle count');
        title('Histogram to show particle speed');
        
        subplot(3,1,2);
        plot(state(1:ecount,1)./1e-9,state(1:ecount,2)./1e-9,'o')
        xlim([0 ConductorL/1e-9])
        ylim([0 ConductorW/1e-9])
        xlabel('Electrons x position (nm)')
        ylabel('Electrons y position (nm)')
        title('Electrons with collision and Mean Free Path')
        
        subplot(3,1,3);
        plot(detaT*(0:i-1),temp(1:i));
        xlabel('time(s)');
        ylabel('Temperature (K)');
        title('Temperature of semiconductor over time');
        
     end
    end
    
    %iterate over the visible particles and show thier
    %paths
    figure(4)
    hold on;
    xlim([0 ConductorL/1e-9]);
    ylim([0 ConductorW/1e-9]);
    xlabel('Electrons x position (nm)');
    ylabel('Electrons y position (nm)');
    title('Trajectories of Electrons with after collisions and MFP');
    for i = 1: ecount
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '-');
    end
    
    Vx=(vth/sqrt(2))*rand(ecount,1);
    Vy=(vth/sqrt(2))*rand(ecount,1);
    Vdis=sqrt(Vx.^2+Vy.^2);
    Vavg = mean(Vdis);
    MFP = Vavg*tmn;