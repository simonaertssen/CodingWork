
%% Introduction
% This code was developed by Simon Aertssen (s181603) and contains the
% numerical solver LMSolver with the framework needed to reproduce results
% from the report. The script runs through different combinations of a and
% lambda, calculations pause in between to show the results. All functions 
% including the LMsolver can be found at the end of the script.

set(0, 'DefaultLineLineWidth', 1.5);
set(0,'DefaultAxesXGrid','on')
set(0,'DefaultAxesYGrid','on')

%% Data and Setup
t0 = 0; tend = 10;
tspan = [t0,tend];
f1 = figure('Renderer', 'painters', 'Position', [300 300 900 600]);
% Use the following figure when printing results:
%f1 = figure('units','normalized','outerposition',[0 0 1 1])

%% Loop: reproduce the results
% Produce a figure for every combination of a and lambda for all three
% step sizes, a total of 9 figures with 6 subplots. 

% Use plotnmr to save plots correctly
plotnmr = 1;
for lambda = [0,1,20]
    param.lambda = lambda;
    for a = [1,2,10]
        param.a = a;
        y0 = param.a; 
        % Use index to display plots correctly
        index = 1;
        for h = [0.1, 0.05, 0.025]
            % Find number of points from given stepsize h
            n = (tspan(end)-tspan(1))/h;


%% Numerical solution:
[tout, yout] = LMsolver(@(t,y,param) func(t,y,param), tspan, n, y0, param);

%% Analytical solution:
% Run the script with the desired analytical solution:
%exact = (tout.^4 + 2*tout.^2 + 1);
%exact = (tout.^4 + 2*tout.^2 + a);

if lambda == 0
    exact = (tout.^4 + 2*sqrt(a)*tout.^2 + a);
else
    exact = (a-1)*exp(-lambda.*tout) + (tout.^4 + 2*tout.^2 + 1);
end

%% Validation with ode45:
% Use ode45 on the same timespan as LMsolver to validate results
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[t,y] = ode45(@(t,y) func(t,y,param), linspace(t0,tend,n+1), y0, options);
 
%% Plots
subplot(2,3,index)
hold on
% LMSolver solution:
plot(tout,yout,'r') 
% Exact solution: (display as 10 points)
s = round(length(tout)/10);
plot(tout(1:s:end), exact(1:s:end),'ob','Markersize', 7,'LineWidth',1)
% ode45 solution: (display as 10 points)
plot(t(1:s:end),y(1:s:end), '.k', 'MarkerSize', 10,'LineWidth',1) 
xlabel('Time [t]'); ylabel('y(t)'); 
title(['a = ',num2str(param.a),', \lambda = ',num2str(param.lambda),' and h = ',num2str(h)])
legend('Numerical solution', 'Exact solution', 'ode45', 'Location','NorthWest')
hold off

subplot(2,3,index+3)
hold on
E = exact - yout;
plot(tout, E, 'r')
% Validate difference between LMsolver and ode45:
% plot(tout, yout - y, '--k')
xlabel('Time [t]'); ylabel('E = y(t_n) - Y_n'); 
title('Global error')
hold off

index = index + 1;
        end
        
        % Print results if desired:
        %print(['Images/Exact2/Figure',num2str(plotnmr),'.png'], '-dpng', '-r200')
        plotnmr = plotnmr + 1;
        
        % Pause to see the results
        pause
        clf(f1)
    end
end


close all

%% Functions 
function dydt = func(t,y,param)
    % Contains the description of the differental equation
    dydt = 4*t*sqrt(y) - param.lambda*(y - (1 + t^2)^2);
    %f = 4*t^3 + 4*t;
end

function [tout,yout] = LMsolver( func, tspan, n, y0, param)
    % This function integrates a given function over tspan using a fourth order
    % multistep predictor-corrector method.
    
    %% Part 0: declaring lengths and allocating space
    h = (tspan(end)-tspan(1))/n;
    tout = linspace(tspan(1), tspan(end), n+1)';
    yout = zeros(1,length(tspan))';
    yout(1) = y0;
    % Declare tolerance: Matlab has precision of 16 digits
    tol = sqrt(10^(-15));
    
    
    %% Part 1: getting started
    % Integrate first three points with an RK scheme:
    for i = 1:3
        K1 = h*func(tout(i), yout(i), param);
        K2 = h*func(tout(i) + h/2, yout(i) + K1/2, param);
        K3 = h*func(tout(i) + h/2, yout(i) + K2/2, param);
        K4 = h*func(tout(i) + h  , yout(i) + K3,   param);
        yout(i+1) = yout(i) + 1/6*(K1 + 2*K2 + 2*K3 + K4);
    end
    
    %% Part 2: multistep method
    % Continue working with the Adams-Bashforth and Adams-Moulton methods
    for i = 4:n
        % Calculate ytemp with explicit method
        ytemp = yout(i) + h/24*( 55*func(tout(i), yout(i),param) - 59*func(tout(i-1), yout(i-1),param) ...
            + 37*func(tout(i-2), yout(i-2),param) - 9*func(tout(i-3), yout(i-3),param) );
        
        % Iterate over implicit scheme: save last two values for y(i+1) to
        % calculate tolerance and thend discard of results.
        iter = [0, ytemp];
        while abs(iter(end) - iter(1)) >= tol
            iter(1) = iter(end);
            iter(end) = yout(i) + h/24*( 9*func(tout(i+1), ytemp,param) + 19*func(tout(i), yout(i),param)...
                - 5*func(tout(i-1), yout(i-1),param) + 1*func(tout(i-2), yout(i-2),param));
        end
        yout(i+1) = iter(end);
    end
end

