% David Talson
% ELEC4700 Assignment 3
clear all
clc


nx=200;
ny=100;
CurrDen = [];
G=sparse(nx*ny,nx*ny);
B=zeros(nx*ny,1);
x2=180;
x1=20;
y2=60;
y1=40;
condRegion=1;
condBottle=0.01;
cond=1.*ones(nx,ny);

for i = 1:nx
    for j = 1:ny
        if(((i>=x1)&&(i<=x2)&&(j<=y1))||((i>=x1)&&(i<=x2)&&(j>=y2)))
            cond(i,j) = condBottle;
        end
    end
end

for i=1:nx
    for j=1:ny
        
        n = j + (i-1)*ny;
        nxm = (i-2)*ny + j;
        nxp = i*ny + j;
        nym = (i-1)*ny + j-1;
        nyp = (i-1)*ny + j+1;
        
        if i==1
            B(n,1)=1;
            G(n,n)=1;
        elseif i==nx
            B(n,1)=0;
            G(n,n)=1;
        elseif j==1
            B(n,1)=0;
            
            condxm = (cond(i,j) + cond(i-1,j))/2;
            condxp = (cond(i,j) + cond(i+1,j))/2;
            condyp = (cond(i,j) + cond(i,j+1))/2;
            
            G(n,n) = -(condxm+condxp+condyp);
            G(n,nxm) = condxm;
            G(n,nxp) = condxp;
            G(n,nyp) = condyp;
            
        elseif j==ny
            B(n,1)=0;
            
            condxm = (cond(i,j) + cond(i-1,j))/2;
            condxp = (cond(i,j) + cond(i+1,j))/2;
            condym = (cond(i,j) + cond(i,j-1))/2;
            
            G(n,n) = -(condxm+condxp+condym);
            G(n,nxm) = condxm;
            G(n,nxp) = condxp;
            G(n,nym) = condym;
        else
            B(n,1) = 0;
            condxm = (cond(i,j) + cond(i-1,j))/2;
            condxp = (cond(i,j) + cond(i+1,j))/2;
            condyp = (cond(i,j) + cond(i,j+1))/2;
            condym = (cond(i,j) + cond(i,j-1))/2;

            G(n,n) = -(condxm+condxp+condyp+condym);
            G(n,nxm) = condxm;
            G(n,nxp) = condxp;
            G(n,nym) = condym;
            G(n,nyp) = condyp;
            
        end
    end  
end
V=G\B;
for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        map(i,j) = V(n,1);
    end
end

[Ex,Ey] = gradient(map',1,1);

figure(5)
surf(cond)
colorbar
grid on
title('Conductance Plot')
xlabel('Region Length')
ylabel('Region Width')

figure(6)
surf(map)
colorbar
grid on 
title('Voltage Plot of Region')
xlabel('Region Length')
ylabel('Region Width')

figure(7) 
quiver(Ex,Ey)
title('Electric Field of Region')
xlabel('Region Width')
ylabel('Region Length')
grid on 

figure(8)
quiver(cond'.*(-Ex),cond'.*(-Ey))
title('Current Density of Region')
xlabel('Region Width')
ylabel('Region Length')
grid on 

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

% Check if there are electrons in the box
outbox = 1;
for chk = 1:numElectrons
        
         while (electronXpos(chk) >= 0.8e-7 & electronXpos(chk) <= 1.2e-7) & ((electronYpos(chk) >= 0 & electronYpos(chk) <= 0.4e-7)|(electronYpos(chk) >= 0.6e-7 & electronYpos(chk) <= 1e-7))
              electronXpos(chk) = rand(1, 1).* regionWidth;
              electronYpos(chk) = rand(1, 1).* regionLength;
         end

end

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
electronAccX =0;
electronAccY =0;
xposNew=electronXpos;
yposNew=electronYpos;
for i=1:1000
      for j = 1:numElectrons 
          if probScat > rand
                angleNew = rand(1).*2*pi;
                electronXvel(1,j) = random(MBdist,1).*cos(angleNew);
                electronYvel(1,j) = random(MBdist,1).*sin(angleNew);
          end
          
          
             if electronYpos(1,j) + electronYvel(1,j).*deltaT >=1e-7|| electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0||((electronYpos(1,j) + electronYvel(1,j).*deltaT >= 0.6e-7||electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0.4e-7) & (electronXpos(1,j) + electronXvel(1,j)*deltaT>=0.8e-7 & electronXpos(1,j) + electronXvel(1,j)*deltaT<=1.2e-7))
                electronYvel(1,j) = electronYvel(1,j)*-1;          
             end
             
             
             yposNew(i,j) = electronYpos(1,j) + electronYvel(1,j).*deltaT;
             if electronXpos(1,j) + electronXvel(1,j)*deltaT >= regionWidth || electronXpos(1,j) + electronXvel(1,j)*deltaT <= 0||((electronYpos(1,j) + electronYvel(1,j).*deltaT >= 0.6e-7||electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0.4e-7) & (electronXpos(1,j) + electronXvel(1,j)*deltaT>=0.8e-7 & electronXpos(1,j) + electronXvel(1,j)*deltaT<=1.2e-7))
                if electronXpos(1,j) + electronXvel(1,j)*deltaT >= regionWidth ||  electronXpos(1,j) + electronXvel(1,j)*deltaT >= 0.8e-7
                    xposNew(i,j) = 0;
                else
                    xposNew(i,j) = regionWidth;
                end
             else
                 electronXvel(1,j) = electronXvel(1,j) + electronAccX*(deltaT);
                 xposNew(i,j) = electronXpos(1,j) + electronXvel(1,j).*deltaT + 0.5*electronAccX*(deltaT^2);
             end
              t1 = sqrt(electronXvel(1,j).^2 + electronYvel(1,j).^2);
              temperature(i,j) = ((mean(t1)^2)*C.m_eff)/(2*C.kb);
              tempx(i,j) = electronXpos(1,j);
              tempy(i,j) = electronYpos(1,j);
              
              xii = ceil((10e8)*xposNew(i,j));
              if xii==0
                  xii=1;
              end
              if xii>200
                  xii=200;
              end
         
              yii=ceil((10e8)*yposNew(i,j));
              if yii==0
                  yii=1;
              end
              if yii>100
                  yii=100;
              end
              EX=(10^9)*Ex(yii,xii);
              EY=(10^9)*Ey(yii,xii);
              electronAccX = ((EX)*C.q_0)/C.m_0;
              electronAccY = ((EY)*C.q_0)/C.m_0;      
               t2(i) = (mean(t1)*C.q_0*(10e15)*regionLength*regionWidth); 
        end
    
    electronXpos = xposNew(i,:);
    electronYpos = yposNew(i,:);
end

figure(9)
plot(xposNew,yposNew,'.','MarkerSize',8);
hold on
xlim([-0.1e-7 2.1e-7]);
ylim([-0.1e-7 1.1e-7]);
grid on
xlabel('Electrons x position (m)')
ylabel('Electrons y position (m)')
title('A plot of path of 1000 electrons in random motion with scaterring probability')
hold off

electron_density = [xposNew(1000,:)',yposNew(1000,:)'];

figure(10)
hist3(electron_density,'CDataMode','auto','FaceColor','interp')
xlabel('Electron x position(m)')
ylabel('Electron y position(m)')
title('Electron Density Map of final positions of electron')
colorbar

figure(11)
plot(t2,'Linewidth',2)
title('Current over time in x direction')
xlabel('Simulation time')
ylabel('Current')
grid on

