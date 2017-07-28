clear variables;
close all
clc

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

dx = 1e-6 %cell width in meter
   %Armin: from the text 
N = 100 %size of the simulation domain in cells
   %Armin: from the text 
dt = 1e-3 %time step in seconds
   %Armin: from the text 
Nt = 1000 %number of timesteps
   %Armin: from the text 
Pi = 3.1415; %acos(-1)
% C0 = ..; %initial concentration of phase 0 in at.%
% C1 = ..; %initial concentration of phase 1 in at.%
epsilon = 1; %constant, forming the coefficient of the non-local energy contribution
   %Armin: from the text 
W = 1e12;
   %Armin: from the text 
M = 1e-10; %interface mobility in m^4/(Js)
   %Armin: from the text 
intfWidth = 1e-7; %interface width in cells
   %Armin: I chooze the number!!
intfPosition = 1e-6/2; %index of the cell where the interface is located
   %Armin: middle of the cell!!
% D = ..; %relaxation coefficient of the diffusion equation

%discretised space containing N cells, using dx as space discretisation
x = 0:dx:(N-1)*dx;

%setting the phase field variable to initial values
phi(1,1:(N/10)) = 0; 
   % Armin: first half of the cell
phi(1,(N/10)+1:N) = 1; %where should phi be equal to 1?
   % Armin: 2nd half of the cell
%setting initial concentration using the initial phase composition values
 %c(1, 1:N/2) = ..;

%time loop
for step = 2:Nt
    %using previous state of the phase field variable as starting condition
    %in the current step
    phi(step, 1:N) = phi(step-1, 1:N);
%     c(step, 1:N) = c(step-1, 1:N);
    %preparing the storage structure for laplace(phi) and laplace(c)
    laplacePhi = x*0;
%     laplaceMu = x*0;
    
    %setting the driving force
    dG = G1const(1) - G0const(1);
       %Armin: from the text (G1-G0)
    
    %setting the chemical potential
%     mu = ..;
    
    %calculating the second spatial derivative of phase field and chemical
    %potential, called laplacePhi and laplaceMu in the code
    for i = 2:N-1
        laplacePhi(i) = (phi(step,i+1) - 2*phi(step,i) + phi(step,i-1))/(dx*dx);
%         laplaceMu(i) = (mu(i+1) - 2*mu(i) + mu(i-1))/(dx*dx);
    end

    %calculating phase field and concentration increments using the given
    %differential equations
    
    div2phi= zeros(1,N);
       % Armin: divergence^2 phi, numerical method from the text !! check again 
    for i= 2 : N-1
      div2phi(1,i)= (1/(dx^2))*(phi(step, i-1)-2*phi(step, i)+phi(step, i+1));
    end   
       % Armin: from the text !! check again

    phidot(step, 1:N) = M*( epsilon^2*div2phi - gPrime(phi(step,:))*W ...
                        + g(phi(step,:))*dG*36/(20*dx));
             
             
             
        % Armin: calculation of phidor from the text check again!!!

%     cdot = ..;
        
    %integrating the phase field and concentration in time from t=(step-1)*dt
    %to t=step*dt
    phi(step, 1:N) = phi(step-1,1:N) + dt * phidot(step-1,1:N);
        %Armin: from the text
     
%     c(step, 1:N) = ..;
    
    %setting boundary conditions phase field and for concentration
    phi(step, 1) = phi(step, 2);
    phi(step, N) = phi(step, N-1);
%     c(step, 1) = c(step, 2);
%     c(step, N) = c(step, N-1);
    
    %plotting
    if(mod(step,Nt/10) == 0 || step == 2)
        plot(phasePlot, x, phi(step, 1:N));
        plot(phaseDotPlot, x, phidot);
        plot(dGPlot, x, dG);
        
%         plot(concPlot, x, c(step, 1:N));
%         plot(conc0Plot, x, c0(step, 1:N));
%         plot(conc1Plot, x, c1(step, 1:N));
    end;
end

%++++++++++++++++++++++++++++++++auxilary functions++++++++++++++++++++++++++++%
%double well function
function [gPhi] = g(phi)
  gPhi = ( (phi(1,:).* phi(1,:)) ) .* ( (1-phi(1,:)) .* (1-phi(1,:)) ) ;
end

%derivative of double well function with respect to phi
function [gPrimePhi] = gPrime(phi)
   gPrimePhi = 2 * phi(1,:) .* ( 2 * ( phi(1,:) .* phi(1,:) ) - 3 * ...
               ( phi(1,:) ) + 1 );
end

%constant phase gibbs energy of phase 1
function [G0c] = G0const(c)
   G0c = 0;
end

%constant phase gibbs energy of phase 0
function [G1c] = G1const(c)
   G1c = -1e6;
end

% %phase gibbs energy of phase 0
% function [G0c] = G0(c)
%  %%%%   G0c = ..;
% end
% 
% %phase gibbs energy of phase 1
% function [G1c] = G1(c)
%     G1c = ..;
% end
% 
% %chemical potential of phase 0
% function [MUc] = mu0(c)
%     MUc = ..;
% end
% 
% %phase gibbs energy of phase 1
% function [MUc] = mu1(c)
%     MUc = ..;
% end