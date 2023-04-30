clc;clear;close all;

profile on
t_range = linspace(0,32.9,258);
df_max = 258/32.9/2;
df_range = linspace(-df_max,df_max,258);
delay_range = linspace(0,128,2048);
fai_range = linspace(-8,8,2048);
load('.\NOF1\mat\NOF1_001.mat');
%% calculate h
h_abs = mapminmax(abs(h),0,1);
h_log = 20*log10(h_abs);

%% calculate S
size_h = size(h);
s = zeros(size(h));
for i = 1:size_h(2)
    s(:,i) = abs(fftshift(fft(h(:,i))));
end
s_map = s/max(s,[],'all');
% s_map = mapminmax(abs(s),0,1);
s_log = 20*log10(s_map);

%% calculate H
H = zeros(size_h);
for i = 1:size_h(1)
   H(i,:) = abs(fftshift(fft(h(i,:))));
end
% H_map = H/max(H,[],'all')+0.0000001;
H_map = mapminmax(abs(H),0,1);
H_log = 20*log10(H_map);
%% calculate B 
B = zeros(size_h);
for i = 1:size_h(1)
    B(i,:) = abs(fftshift(fft(s(i,:))));
end
% 20*log10(mapminmax(abs(B),0,1));
B_map = mapminmax(abs(B),0,1);
B_log = 20*log10(B_map);

%% calculate Rs(0,phi)
Doppler_sp = zeros(1,size_h(1));
for i = 1:size_h(1)
   Doppler_sp(i) = sum(s(i,:).^2); 
end
Doppler_sp_abs = abs(Doppler_sp);
Doppler_sp_map = Doppler_sp_abs/max(Doppler_sp);
Doppler_sp_log = 10*log10(Doppler_sp_map);

%% calculate Rh(tau,0)
p = zeros(1,2048);
for i = 1:2048
    p(i) = (sum(s(:,i).^2));
end
p_map = p/max(p);
p_log = 10*log10(p_map);

%% DRAW Transfer function
figure(1)
tx1 = suptitle("NOF\_001's System Functions");
set(tx1,'position',get(tx1,'position')+[0 0.02 0]);
subplot(231)
load clown
image(delay_range,t_range,h_log,'CDataMapping','scaled');
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
xticks(0:16:128);
title("$10log_{10}(|\hat{h}(\tau,t)|^2)$",'interpreter','latex');
xlabel("Time Delay $\tau$/ms",'interpreter','latex');
ylabel("Time t/ms",'interpreter','latex');

subplot(232)
load clown
image(fai_range,t_range,H_log,'CDataMapping','scaled')
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
xticks(-8:2:8);
title("$10log_{10}(|\hat{H}(f,t)|^2)$",'interpreter','latex');
xlabel("Frequency f(kHz)",'interpreter','latex');
ylabel('Time t/ms','interpreter','latex');

subplot(234)
load clown
image(delay_range,df_range,s_log,'CDataMapping','scaled');
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
xticks(0:16:128);
yticks(-4:1:4);
title("$10log_{10}(|\hat{s}(\tau,\varphi)|^2)$",'interpreter','latex');
xlabel("Time delay $\tau$(ms)",'interpreter','latex');
ylabel("Frequency Shift $\varphi$ (Hz)",'interpreter','latex');

subplot(235)
load clown
image(fai_range,df_range,B_log,'CDataMapping','scaled')
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
xticks(-8:2:8);
yticks(-4:1:4);
title("$10log_{10}(|\hat{B}(f,\varphi)|^2)$",'interpreter','latex');
xlabel("Frequency f(kHz)",'interpreter','latex');
ylabel('Frequency Shift $\varphi$ (Hz)','interpreter','latex');

subplot(233)
plot(Doppler_sp_log,df_range,'b')
xlim([-40,0])
title("$10log_{10}(|\hat{s}(0,\varphi)|^2)$",'interpreter','latex');
xlabel("Power Density (dB)",'interpreter','latex');
ylabel("Frequecy Shift $\varphi$ (Hz)",'interpreter','latex');

subplot(236)
plot(delay_range,p_log,'b');
ylim([-50,0])
title("Power Delay Profile (PDP): $10log_{10}(\int s(\tau,\varphi)d\varphi)$",'interpreter','latex');
xlabel("Time Delay $\tau$(ms)",'interpreter','latex');
ylabel("Power Density(dB)",'interpreter','latex');
drawnow;

%% calculate CW signal in time-varible channel
dt = 128/1000/2048;
fs = 1/dt;
signal_width = 0.05;
signal_t = 0:dt:signal_width;
f_signal_CW = 1000;
rec_t = 0:dt:(signal_width+delay_range(end)/1000)-dt;
corr_t = -(signal_width+delay_range(end)/1000)+dt:dt:(signal_width+delay_range(end)/1000)-dt;

CW_signal = exp(-1j*2*pi*f_signal_CW*signal_t);
rec_signal_CW = zeros(size_h(1),length(rec_t));
for i = 1:size_h(1)
   rec_signal_CW(i,:) = conv(h(i,:),CW_signal);
end
rec_signal_CW_envp = 20*log10(abs(rec_signal_CW)/abs(max(rec_signal_CW,[],'all')));

rec_corr_CW = zeros(size_h(1),length(corr_t));
for i = 1:size_h(1)
   rec_corr_CW(i,:) = abs(xcorr(real(rec_signal_CW(i,:)),real(CW_signal)));
end

%% Draw signal in time-varible channel
figure(2)
subplot(211)
plot(signal_t,real(CW_signal));
subplot(223)
waterfall(rec_t,t_range,rec_signal_CW_envp);
shading interp
subplot(224)
waterfall(corr_t,t_range,rec_corr_CW);
drawnow;

%% calculate LFM signal in time-varible channel
dt = 128/1000/2048;
fs = 1/dt;
signal_width = 0.05;
signal_t = 0:dt:signal_width;
f_signal_LFM_start = 1000;
f_signal_LFM_end   = 1500;
k = (f_signal_LFM_end-f_signal_LFM_start)/signal_width;
rec_t = 0:dt:(signal_width+delay_range(end)/1000)-dt;
corr_t = -(signal_width+delay_range(end)/1000)+dt:dt:(signal_width+delay_range(end)/1000)-dt;

LFM_signal = exp(-1j*2*pi*(f_signal_LFM_start+k*signal_t).*signal_t);
rec_signal_LFM = zeros(size_h(1),length(rec_t));
for i = 1:size_h(1)
   rec_signal_LFM(i,:) = conv(h(i,:),LFM_signal);
end
rec_signal_LFM_envp = 20*log10(abs(rec_signal_LFM)/abs(max(rec_signal_LFM,[],'all')));

rec_corr_LFM = zeros(size_h(1),length(corr_t));
for i = 1:size_h(1)
   rec_corr_LFM(i,:) = abs(xcorr(real(rec_signal_LFM(i,:)),real(LFM_signal)));
end

%% Draw signal in time-varible channel
figure(3)
subplot(211)
plot(signal_t,real(LFM_signal));
subplot(223)
waterfall(rec_t,t_range,rec_signal_LFM_envp);
subplot(224)
waterfall(corr_t,t_range,rec_corr_LFM);
drawnow;
profile viewer
