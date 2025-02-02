clear variables; close all; clc
%% setting up the plots for phase field
        figure
        
        phasePlot = subplot(4,2,1);
        hold on;
        title('phase field');
        ylim([0 1.1])
        
        phaseDotPlot = subplot(4,2,3);
        hold on;
        title('phase field increment rate');
        ylim([-10 10])
        
        dGPlot = subplot(4,2,5);
        hold on;
        title('dG');
        
        G0Plot = subplot(4,2,7);
        hold on;
        title('G0');
                
        %setting up the plots for concentration
        
        concPlot = subplot(4,2,2);
        hold on;
        title('concentration');
        
        conc0Plot = subplot(4,2,4);
        hold on;
        title('phase concentration 0');
        
        conc1Plot = subplot(4,2,6);
        hold on;
        title('phase concentration 1');
        
        G1Plot = subplot(4,2,8);
        hold on;
        title('G1');
%%
dx = 1e-6; %cell width in meter
N = 100; %size of the simulation domain in cells
dt = 1e-3; %time step in seconds
Nt = 20000; %number of timesteps
Pi = 3.1415; %acos(-1)
C0 = 13; %initial concentration of phase 0 in at.%
C1 = 20; %initial concentration of phase 1 in at.%
epsilon = 1; %constant, forming the coefficient of the non-local energy 
             %contribution
W = 1e12;
M = 1e-10; %interface mobility in m^4/(Js)
intfWidth = 20*dx; %interface width in cells
intfPosition = 0.9; %index of the cell where the interface is located
D = 1e-15; %relaxation coefficient of the diffusion equation

%discretised space containing N cells, using dx as space discretisation
x = 0:dx:(N-1)*dx;

%setting the phase field variable to initial values
phi(1,1:(N*intfPosition)) = 0; 
phi(1,(N*intfPosition)+1:N) = 1; %where should phi be equal to 1?

%setting initial concentration using the initial phase composition values
c(1, 1:N) = phi * C1 + (1-phi) * C0;

%% time loop
for step = 2:Nt
    %using previous state of the phase field variable as starting condition
    %in the current step
    phi(step, 1:N) = phi(step-1, 1:N);
    c(step, 1:N)   = c(step-1, 1:N);
    %preparing the storage structure for laplace(phi) and laplace(c)
    laplacePhi = x*0;
    laplaceMu  = x*0;
    
    %setting the driving force
    dG =  G0(c(step, 1:N)) -  G1(c(step, 1:N));    

    %setting the chemical potential
    mu = phi(step, 1:N) .* mu1(c(step, 1:N)) + (1-phi(step, 1:N)) .*...
         mu0(c(step, 1:N));

    %calculating the second spatial derivative of phase field and chemical
    %potential, called laplacePhi and laplaceMu in the code
    for i = 2:N-1
      laplacePhi(i) = (phi(step,i+1) - 2*phi(step,i) + phi(step,i-1))/...
                      (dx*dx);
      laplaceMu(i) = (mu(i+1) - 2*mu(i) + mu(i-1))/(dx*dx);
    end

    %calculating phase field and concentration increments using the given
    %differential equations
    phidot = M * (epsilon^2 * laplacePhi - gPrime(phi(step,:))* W + ...
             g(phi(step,:)) .* dG * 36 / intfWidth);
    cdot = D * laplaceMu;
    
    %integrating the phase field and concentration in time from 
    %t=(step-1)*dt to t=step*dt
    phi(step, 1:N) = phi(step-1,1:N) + dt * phidot;
    c(step, 1:N)   = c(step-1,1:N)   + dt * cdot ;
    
    %setting boundary conditions phase field and for concentration
    phi(step, 1) = phi(step, 2);
    phi(step, N) = phi(step, N-1);
    c(step, 1) = c(step, 2);
    c(step, N) = c(step, N-1);
    
    %plotting
    if(mod(step,Nt/10) == 0 || step == 2)
        plot(phasePlot, x, phi(step, 1:N));
        plot(phaseDotPlot, x, phidot);
        plot(dGPlot, x, dG);
        plot(G0Plot, x, G0(c(step, 1:N)));
        plot(G1Plot, x, G1(c(step, 1:N)));
        plot(concPlot, x, c(step, 1:N));
     end;
end

%+++++++++++++++++++++++++++auxilary functions++++++++++++++++++++++++++++%
%double well function
function [gPhi] = g(phi)
  gPhi = ( (phi(1,:).* phi(1,:)) ) .* ( (1-phi(1,:)) .* (1-phi(1,:)) ) ;
end

%derivative of double well function with respect to phi
function [gPrimePhi] = gPrime(phi)
  gPrimePhi = 2 * phi(1,:) .* ( 2 * ( phi(1,:) .* phi(1,:) ) - ...
              3 * ( phi(1,:) ) + 1 );
end

 function [G0c] = G0(x)
   G0c = ( (13e4*(x-10).*(x-10))/2 + 1e6);
 end
 
 %phase gibbs energy of phase 1
 function [G1c] = G1(y)
   G1c = ( (5e4*(y-20).*(y-20))/2 + 1e5);
 end
 
%chemical potential of phase 0
 function [MUc] = mu0(c)
     MUc = 13e4 * (c-10);
 end
 
%phase gibbs energy of phase 1
 function [MUc] = mu1(c)
     MUc = 5e4 * (c-20);
 end
 
