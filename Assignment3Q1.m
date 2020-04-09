% David Talson
% ELEC4700 Assignment 3
clear all
clc

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
C.g = 9.80665;                      % metres (32.1740 ft) per s²
C.m_eff = 0.26*C.m_0;               % Effective mass of electrons 

regionWidth = 2e-7;                 % Nominal width of the region
regionLength = 1e-7;                % Nominal length of the region
T = 300;                            % Assumed temperature

vth = sqrt((2*C.kb*T)/C.m_eff); % Thermal velocity as mean of magnitude of velocity
tmn = 2e-13;                        % Mean time between collisions
freepath = vth*2e-13;

numElectrons = 1000;
electronXpos = rand(1,numElectrons).*regionWidth;
electronYpos = rand(1,numElectrons).*regionLength;
angle = rand(1,numElectrons).*2*pi;

sig = sqrt(C.kb*T/C.m_eff)/4;
MBdist = makedist('Normal',vth,sig);
electronVel = random(MBdist,1,numElectrons);

angle = rand(1,numElectrons).*2*pi;
electronXvel = electronVel.*cos(angle);
electronYvel = electronVel.*sin(angle);
deltaT = 1e-9/vth;

probScat = 1 - exp(-deltaT/tmn);
electronVel = sqrt(electronXvel.^2 + electronYvel.^2);
electronAcc =-((1.5/regionWidth)*C.q_0)/C.m_0;
count = 0;

for i=1:100
      for j = 1:numElectrons 
            if probScat > rand
                count = count + 1;
                angleNew = rand(1).*2*pi;
                electronXvel(1,j) = random(MBdist,1).*cos(angleNew);
                electronYvel(1,j) = random(MBdist,1).*sin(angleNew);  
            end
            
            if electronYpos(1,j) + electronYvel(1,j).*deltaT  >= regionLength || electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0
               electronYvel(1,j) = electronYvel(1,j)*-1;  
            end
            yposNew(i,j) = electronYpos(1,j) + electronYvel(1,j).*deltaT ;
            
             electronXvel(1,j) = electronXvel(1,j) + electronAcc*(deltaT);
             xposNew(i,j) = electronXpos(1,j) + electronXvel(1,j).*deltaT + 0.5*electronAcc*(deltaT^2);
              
            if xposNew(i,j) >= regionWidth || xposNew(i,j) <= 0
                if xposNew(i,j) >= regionWidth
                    xposNew(i,j) = 0;
                else
                    xposNew(i,j) = regionWidth;
                end
            end
               
           % To calculate and plot the semiconductor temperature
           t1 = sqrt(electronXvel(1,:).^2 + electronYvel(1,:).^2);
           temperature(i,j) = ((mean(t1)^2)*C.m_eff)/(2*C.kb);
           tempx(i,j) = electronXpos(1,j);
           tempy(i,j) = electronYpos(1,j);
    
           %t2(i,j) = (mean(t1)*C.q_0*(10e15)*regionLength*regionWidth);  
           t2(i) = (mean(t1)*C.q_0*(10e15)*regionLength*regionWidth); 
           
        end 
   
    electronXpos = xposNew(i,:);
    electronYpos = yposNew(i,:);
end
mfp = (1000/count)*deltaT*mean(electronVel);
meantime = deltaT*(1000/count);

figure (1)
plot(xposNew,yposNew,'.','Markersize',8)
xlim([-0.1e-7 2.1e-7])
ylim([-0.1e-7 1.1e-7])
grid on
xlabel('Electrons x position (m)')
ylabel('Electrons y position (m)')
title('A plot of path of 50 electrons in random motion with scaterring probability')

figure(2)
plot(t2,'Linewidth',2)
title('Current over time in x direction')
xlabel('Simulation time')
ylabel('Current')
grid on

electron_density = [xposNew(100,:)',yposNew(100,:)'];

figure(3)
hist3(electron_density,'CDataMode','auto','FaceColor','interp')
xlabel('Electron x position(m)')
ylabel('Electron y position(m)')
title('Electron Density Map of final positions of electron')
colorbar


figure(4)
[X,Y] = meshgrid(0e-7:0.1e-7:2e-7,0e-7:0.1e-7:1e-7);
tempNew = griddata(tempx,tempy,temperature,X,Y);
surf(X,Y,tempNew)
colorbar
xlabel('Region Width')
ylabel('Region Length')
title('Plot of temperature map')



