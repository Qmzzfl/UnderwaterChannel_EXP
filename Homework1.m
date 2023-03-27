clear;clc;close all;
%%  parameters setting
%%%%%Determining signal parameters%%%%%
f = 100;
phi = 0.1;
E0 = 10;

%%%%%Simulation parameters%%%%%%%%%%%%
t_end = 1;
fs = 4000;
dt = 1/fs;
t = 0:dt:t_end;

%%%%%%%%%Noise parameters%%%%%%%%%%%%%
SNR = -3;
mag_noise = E0/(10^(SNR/20));
power_noise = mag_noise^2;
%noise generator
u = normrnd(0,sqrt(power_noise/2),1,length(t));
v = normrnd(0,sqrt(power_noise/2),1,length(t));
determine_P = (E0).*cos(2*pi*f*t+phi);
nondetermine_P = u.*cos(2*pi*f*t+phi)-v.*sin(2*pi*f*t+phi);
%% %%%%%%%%Signal Generate%%%%%%%%%%%%%%%%
P = determine_P+nondetermine_P;
%Calculate recieve power and its power spectrum
power_rec = P.^2;
xcor_rec = xcorr(power_rec,'unbiased');
rec_power_sp = fft(xcor_rec,length(t));
rec_sp = abs(rec_power_sp);
%Calculate determine signal and its power spectrum
power_determine = determine_P.^2;
xcor_deter = xcorr(power_determine,'unbiased');
determine_power_sp = fft(xcor_deter,length(t));
determine_sp = abs(determine_power_sp);
%Calculate nondetermine signal and its power spectrum
power_nondetermine = nondetermine_P.^2;
xcor_nondetermine = xcorr(power_nondetermine,'unbiased');
nondetermine_power_sp = fft(xcor_nondetermine,length(t));
nondetermine_sp = abs(nondetermine_power_sp);

%% Draw pic
%Draw signal
figure(1)
plot(t,P);
title('Received Signal(SNR = -3dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
% Calculate xlabel of spectrum
%the midpoint of spectrum corresponds the half of 
df = fs/(4*((length(rec_sp)-1)/2)); 
xlabel_fft = df:df:fs/4;
%Draw power spectrum
figure(2)
subplot(321)
plot(t,P);
title('Received Signal(SNR = -3dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
subplot(322)
plot(xlabel_fft,rec_sp(2:round(length(rec_sp)/2)));
title('Received Signal Power Spectrum(SNR = -3dB)');
xlabel('Frequency/Hz');
ylabel('Power/Pa^2*Hz^-1');
subplot(323)
plot(t,determine_P);
ylim([-50,50])
title('Determine Component of Received Signal(SNR = -3dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
subplot(324)
plot(xlabel_fft,determine_sp(2:round(length(determine_sp)/2)));
title('Determining Signal Power Spectrum(SNR = -3dB)');
xlabel('Frequency/Hz');
ylabel('Power/Pa^2*Hz^-1');
subplot(325)
plot(t,nondetermine_P);
title('Nondetermine Component of Received Signal(SNR = -3dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
subplot(326)
plot(xlabel_fft,nondetermine_sp(2:round(length(nondetermine_sp)/2)));
title('Nondetermining Signal Power Spectrum(SNR = -3dB)');
xlabel('Frequency/Hz');
ylabel('Power/Pa^2*Hz^-1');

figure(3)
subplot(311)
plot(t,P);
title('Received Signal(SNR = -3dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
ylim([-50,50])
subplot(312)
SNR = 10;
mag_noise = E0/(10^(SNR/20));
power_noise = mag_noise^2;
%noise generator
u = normrnd(0,sqrt(power_noise/2),1,length(t));
v = normrnd(0,sqrt(power_noise/2),1,length(t));
P = (E0+u).*cos(2*pi*f*t+phi)-v.*sin(2*pi*f*t+phi);
%Draw signal
plot(t,P);
title('Received Signal(SNR = 10dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
ylim([-50,50])
subplot(313)
SNR = 40;
mag_noise = E0/(10^(SNR/20));
power_noise = mag_noise^2;
%noise generator
u = normrnd(0,sqrt(power_noise/2),1,length(t));
v = normrnd(0,sqrt(power_noise/2),1,length(t));
P = (E0+u).*cos(2*pi*f*t+phi)-v.*sin(2*pi*f*t+phi);
%Draw signal
plot(t,P);
title('Received Signal(SNR = 40dB)');
xlabel('Time/s');
ylabel('Magnitude/Pa');
ylim([-50,50])