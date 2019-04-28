%% Question 1
%Simon Aertssen
clear all
close all
clc
format short 
%Large amounts of numbers, so we don't want them to be very long upon
%inspection. For convenience, all results were plotted in h5 at the bottom,
%still closed at this time.

%% Data and Setup
% websave('heli_signal.txt','https://tinyurl.com/BVC-heli-noise')
%I prefer to present the given data on rows, to have it take less space. 
T = 0.5; Fs = 2000; dt = 1/Fs; t = 0:dt:(T-dt);
ynoise = load('heli_signal.txt'); N = length(ynoise);

%% 1. Plot of the original signal
h1 = figure(1)
plot(t,ynoise)
grid on; box on;
xlabel('Time t [s]'); ylabel('Pressure [N/m^2]')
title('Original signal')
% close(h1)

%% 2. Plot of the FFT
Nq = N/2; f = (0:Nq-1)/T;
Ynoise = abs(fft(ynoise))/Nq; Ynoise = Ynoise(1:Nq);

h2 = figure(2)
plot(f,Ynoise)
grid on; box on;
xlabel('Frequency f [Hz]'); ylabel('Absolute value')
title('FFT(signal)')
% close(h2)

%% 3. Noise levels
%Upon visual inspection from the last figure,, the noise threshold should
%be at least 10. Since the spikes are very large in relation to the noise,
%we will opt for a threshold of 12. At the same time, we will save the
%spike frequencies in F.
F = [];
for i = 1:Nq
    if Ynoise(i) > 12
        Y(i) = Ynoise(i);
        F(end+1) = Ynoise(i);
    else
        Y(i) = 0;
    end
end

%% 4. Plot of the filtered FFT
h3 = figure(3)
plot(f,Y)
grid on; box on;
xlabel('Frequency f [Hz]'); ylabel('Absolute value')
title('FFT(signal) - filtered')
% close(h3)

%% 5. Plot of the filtered signal
y = real(ifft(Y))*Nq; 

h4 = figure(4)
plot((0:dt*2:0.5-dt*2),y)
grid on; box on;
xlabel('Time t [s]'); ylabel('Pressure [N/m^2]')
title('Filtered signal')
% close(h4)

%% 6. Rotations of the rotors
%The Fourier analysis points out that the main rotor turns at a freq of
%80Hz. For three blades, that equals:
Main_rot_speed = F(2)/3;

%For the tip rotor, at a freq of 50 Hz and two blades, this gives:
Tip_rot_speed  = F(1)/2;

%For each rotation of the main rotor, the tip rotor turns 2.4 times

%% 7. Summation of the results
%There are other ways of plotting multiple graphs in one plot, but in this
%way I can quickly turn one 'off' if needed.
h5 = figure(5)
subplot(2,1,1)
hold on
plot(t,ynoise,'k')
plot((0:dt*2:0.5-dt*2),y, 'r')  
hold off; grid on; box on;
xlabel('Time t [s]'); ylabel('Pressure [N/m^2]')
title('Time domain')
legend('Original signal', 'Filtered signal', 'Location', 'southeast')

subplot(2,1,2)
hold on
plot(f,Ynoise)
plot(f,Y)
hold off; grid on; box on;
xlabel('Frequency f [Hz]'); ylabel('Absolute value')
title('Frequency domain')
legend('Original signal', 'Filtered signal', 'Location', 'southeast')
close(h5)

clc
fprintf('For each rotation of the main rotor, the tip rotor turns %.4f times\n', (Tip_rot_speed/Main_rot_speed))


