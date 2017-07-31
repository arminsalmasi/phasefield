clear variables; clc; close all;

%% Armin
% 1- Please check my comments in 1a, 1b and 2a
% 2- phase one grows from left to right
% 3- It looks like that change in c, c0 and c1 follows the grow
% 4- with Nt=10 (small number of time steps) it is possible to watch the 
%   growing phase 
% 5- it looks like there is an impingment at Nt = 20000, the interface
%   passes the cell limit
%%

%setting up the plots for phase field
figure
phasePlot = subplot(3,2,1);
hold on;
title('phase field');
ylim([0 1.1])
phaseDotPlot = subplot(3,2,3);
hold on;
title('phase field increment rate');
ylim([-10 10])
dGPlot = subplot(3,2,5);
hold on;
title('dG');
% ylim([-5e10 5e10])

%setting up the plots for concentration
concPlot = subplot(3,2,2);
hold on;
% ylim([9 21])
title('concentration');
conc0Plot = subplot(3,2,4);
hold on;
% ylim([9 21])
title('phase concentration 0');
conc1Plot = subplot(3,2,6);
hold on;
% ylim([9 21])
title('phase concentration 1');

dx = 1e-6; %cell width in meter
N = 100; %size of the simulation domain in cells
dt = 1e-3; %time step in seconds
Nt = 20000; %number of timesteps
intfWidth = 20*dx; %interface width in cells
Pi = 3.1415; %acos(-1)
m01 = 0.7; %slope of the equilibrium line of phase 0 in contact with phase 1 in the linear phase diagram
m10 = 1; %slope of the equilibrium line of phase 1 in contact with phase 0 in the linear phase diagram
cS = 25; %intersection concentration of the linearised phase diagram
C0 = 13; %initial concentration of phase 0 in at.%
C1 = 20; %initial concentration of phase 1 in at.%
sigma = 1; %interface energy in J/m�
M = 1e-10; %interface mobility in m^4/(Js)
intfPosition = 0.9; %index of the cell where the interface is located %% # Armin: Where is it used?
D = 1e-15; %relaxation coefficient of the diffusion equation

c0(1:Nt, 1:N) = 0;
c1 = c0;

%discretised space containing N cells, using dx as space discretisation
x = 0:dx:(N-1)*dx;

%setting the phase field variable to initial values
phi(1, 1:N*intfPosition) = 0;
phi(1, N*intfPosition+1:N) = 1; %where should phii be equal to 1?

%setting initial concentration using the initial phase composition values
 c(1, 1:N) = phi * C0 + (1-phi) * C1;
%setting initial phase concentration using the partitioning function
 [c0(1, 1:N), c1(1, 1:N)] = partitioning(phi(1, 1:N), c(1, 1:N), m01, m10, cS);

%time loop
for step = 2:Nt
    %using previous state of the phase field variable as starting condition
    %in the current step
    phi(step, 1:N) = phi(step-1, 1:N);
    c(step, 1:N) = c(step-1, 1:N);
    %preparing the storage structure for laplace(phi) and laplace(c)
    laplacePhi = x*0;
    laplaceMu = x*0; 
    
    %setting the driving force
    dG = abs(G1(c1(step, 1:N)) - G0(c0(step, 1:N)));
    
    %setting the chemical potential
     mu = mu1(c1(step, 1:N))-mu0(c0(step, 1:N));
    
    %calculating the second spatial derivative of phase field and chemical
    %potential, called laplacePhi and laplaceMu in the code
    for i = 2:N-1
        laplacePhi(i) = (phi(step,i+1) - 2*phi(step,i) + phi(step,i-1))/(dx*dx);
        laplaceMu(i) = (mu(i+1) - 2*mu(i) + mu(i-1))/(dx*dx);
    end

    %calculating phase field and concentration increments using the given
    %differential equations
    phidot = M * (sigma * ( laplacePhi - gPrime(phi(step,:)) * 27 / (intfWidth^2) ) + phi(step,:) .* (1 - phi(step,:)) .* dG * 6 / intfWidth);
    cdot = D * laplaceMu;
        
    %integrating the phase field and concentration in time from t=(step-1)*dt
    %to t=step*dt
    phi(step, 1:N) = phi(step-1,1:N) + dt * phidot;
    c(step, 1:N) = c(step-1, 1:N) + dt * cdot;
    
    %setting boundary conditions phase field and for concentration
    phi(step, 1) = phi(step, 2);
    phi(step, N) = phi(step, N-1);
    c(step, 1) = c(step, 2);
    c(step, N) = c(step, N-1);
    
     [c0(step, 1:N), c1(step, 1:N)] = partitioning(phi(step, 1:N), c(step, 1:N), m01, m10, cS);

    %plotting
    if(mod(step,Nt/10) == 0 || step == 2)
        plot(phasePlot, x, phi(step, 1:N));
        plot(phaseDotPlot, x, phidot);
        plot(dGPlot, x, dG);
        
         plot(concPlot, x, c(step, 1:N));
         plot(conc0Plot, x, c0(step, 1:N));
         plot(conc1Plot, x, c1(step, 1:N));
    end;
end

%++++++++++++++++++++++++++++++++auxilary functions++++++++++++++++++++++++++++%

%derivative of double well function with respect to phi
function [gPrimePhi] = gPrime(phi)
    gPrimePhi = 2 * phi(1,:) .* ( 2 * ( phi(1,:) .* phi(1,:) ) - 3 * ( phi(1,:) ) + 1 );
end

%constant phase gibbs energy of phase 1
function [G0c] = G0const(c)
    G0c = c ;
end

%constant phase gibbs energy of phase 0
function [G1c] = G1const(c)
    G1c = c -1e6;
end
% %phase gibbs energy of phase 0
 function [G0c] = G0(c)
     G0c = (13e4 * (c-10) .* (c-10) /2 + 1e6);;
 end
 
 %phase gibbs energy of phase 1
 function [G1c] = G1(c)
     G1c = (5e4 * (c-20) .* (c-20) /2 + 1e5);
 end
% 
% %chemical potential of phase 0
 function [MUc] = mu0(c)
     MUc = 13e4 * (c-10);
 end
% 
 %phase gibbs energy of phase 1
 function [MUc] = mu1(c)
     MUc = 5e4 * (c-20);
 end
% 
%%partitioning function that splits the total concentration into individual phase
%%concentrations
function [pC0, pC1] = partitioning(phi, c, m01, m10, cS)
     pC1 = (c + (1-phi) * cS * ((m01/m10)-1)) ./ (phi + (1-phi) * (m01/m10));
     pC0 = cS + (m01/m10) * (pC1-cS);
end
