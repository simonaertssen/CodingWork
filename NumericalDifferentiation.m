clear all 
close all
clc

% This code contains the implementation of Euler, Euler-Heun and
% Runge-Kutta methods of numerical solutions to differential equations.

for j = 3
%% 1 gegevens en analytische oplossing
l = 0.6; g = 9.81; 
h = 1/(4^j); t = 0:h:6; theta = (pi/18)*cos(sqrt(g/l)*t);


%% 2 Euler
thetaE = [pi/18]; phiE = [0]; errorsE = [0];
for i = 1:length(t)-1
    phiE(i+1) = phiE(i) + h*(-g/l)*(thetaE(i));
    thetaE(i+1) = thetaE(i) + h*(phiE(i));
    
    errorsE(i+1) = abs(thetaE(i) - theta(i));
end

%% 3 Euler-Heun
thetaEH = [pi/18]; phiEH = [0]; errorsEH = [0];
for i = 1:length(t)-1
    phiEH(i+1)  = phiEH(i) + (h/2)*( (-g/l)*(thetaEH(i)) + (-g/l)*(thetaE(i+1)) );
    thetaEH(i+1) = theta(i) + (h/2)*( (phiEH(i)) + phiE(i+1));
    
    errorsEH(i+1) = abs(thetaEH(i) - theta(i));
end

%% 4 Runge - Kutta
thetaRK = [pi/18]; phiRK = [0]; errorsRK = [0];
for i = 1: (length(t) - 1)
    k1ph = h*( (-g/l)*thetaRK(i) );
    k1th = h*( phiRK(i) );
    
    k2ph = h*( (-g/l)*(thetaRK(i) + k1th/2) );
    k2th = h*( phiRK(i) + k1ph/2 );
    
    k3ph = h*( (-g/l)*(thetaRK(i) + k1th/4 + k2th/4) );
    k3th = h*( phiRK(i) + k1ph/4 + k2ph/4 );
    
    k4ph = h*( (-g/l)*(thetaRK(i) - k2th + 2*k3th) );
    k4th = h*( phiRK(i) - k2ph + 2*k3ph );
    
    phiRK(i+1) = phiRK(i)  + (1/6)*( k1ph + 4*k3ph + k4ph );
    thetaRK(i+1) = thetaRK(i) + (1/6)*( k1th + 4*k3th + k4th );
    
    errorsRK(i+1) = abs(thetaRK(i) - theta(i));
end

%% Plots
h1 = figure(1)

subplot(4,1,1)
plot(t, theta, '.r', 'Markersize', 1)
grid on
title(['Analytical solution: h = ' num2str(h)])
xlabel('t')
ylabel('\theta (t)')

subplot(4,1,2)
plot(t, thetaE,'.b', 'Markersize', 1)
grid on
title(['Euler solution: h = ' num2str(h)])
xlabel('t')
ylabel('\theta (t)')

subplot(4,1,3)
plot(t, thetaEH,'.g', 'Markersize', 1)
grid on
title(['Euler - Heun solution: h = ' num2str(h)])
xlabel('t')
ylabel('\theta (t)')

subplot(4,1,4)
plot(t, thetaRK,'.k', 'Markersize', 1)
grid on
title(['Runge - Kutta solution: h = ' num2str(h)])
xlabel('t')
ylabel('\theta (t)')

h2 = figure(2)
hold on
plot(t, theta,  'r', 'Markersize', 1)
plot(t, thetaE, 'b', 'Markersize', 1)
plot(t, thetaEH,'g', 'Markersize', 1)
plot(t, thetaRK,'k', 'Markersize', 1)
hold off
grid on 
title(['Comparison of the different solutions for h = ' num2str(h)])
legend('Analytical', 'Euler', 'Euler-Heun', 'Runge Kutta')
xlabel('t')
ylabel('\theta (t)')

h3 = figure(3)
hold on 
plot(t, errorsE, 'b', 'Markersize', 1)
plot(t, errorsEH,'g', 'Markersize', 1)
plot(t, errorsRK,'k', 'Markersize', 1)
hold off
grid on 
title(['Comparison of the error values for h = ' num2str(h)])
legend('Euler', 'Euler-Heun', 'Runge Kutta')
xlabel('t')
ylabel('error')

print (h2, '-djpeg', ['Comparison of the different solutions for h = ' num2str(h) '.jpeg']);
print (h3, '-djpeg', ['Comparison of the different error values for h = ' num2str(h) '.jpeg']);

end

