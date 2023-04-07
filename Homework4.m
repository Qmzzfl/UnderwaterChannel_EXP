clc;clear;close all;

t_range = linspace(0,32.9,258);
df_max = 258/32.9/2;
df_range = linspace(-df_max,df_max,258);
delay_range = linspace(0,128,2048);
fai_range = linspace(-8,8,2048);
load('.\NOF1\mat\NOF1_001.mat');
%% calculate h
suptitle("NOF\_001的各个时变参数函数");
subplot(231)
h_abs = mapminmax(abs(h),0,1);
h_log = 20*log10(h_abs);

%% calculate S
size_h = size(h);
s = zeros(size(h));
for i = 1:size_h(2)
    s(:,i) = abs(fftshift(fft(h(:,i))));
end
s_map = s/max(s,[],'all');
s_log = 20*log10(s_map);

%% calculate B & RB(0,phi)
B = zeros(size_h);
for i = 1:size_h(1)
    B(i,:) = fftshift(fft(s(i,:)));
end

Rb0_abs = abs(B(:,1025));
Rb0_map = Rb0_abs/max(Rb0_abs);
Rb0_log = 20*log10(Rb0_map);
%% calculate Rh(tau,0)
p = zeros(1,2048);
for i = 1:2048
    p(i) = (sum(s(:,i)))^2;
end

%% DRAW Transfer function
figure(1)
subplot(231)
load clown
image(delay_range,t_range,h_log(1:258,1:4:end),'CDataMapping','scaled');
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
title("$10log_{10}(|\hat{h}(\tau,t)|^2)$",'interpreter','latex');
xlabel("时延{\tau}/ms");
ylabel("绝对时间{t}/ms");

subplot(232)
load clown
image(delay_range,df_range,s_log,'CDataMapping','scaled');
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
title("$10log_{10}(|\hat{s}(\tau,\varphi)|^2)$",'interpreter','latex');
xlabel("时延{\tau}(ms)");
ylabel("多普勒频移{\phi}(Hz)");


subplot(233)
load clown
image(fai_range,df_range,20*log10(mapminmax(abs(B),0,1)),'CDataMapping','scaled')
colormap(slanCM(153,160))
caxis([-40 0])
colorbar
title("$10log_{10}(|\hat{B}(f,\varphi)|^2)$",'interpreter','latex');
xlabel("响应频率{f}(kHz)");
ylabel("多普勒频移{\phi}(Hz)");

subplot(223)
plot(Rb0_log,df_range)
title("$10log_{10}(|\hat{s}(0,\varphi)|^2)$",'interpreter','latex');
xlabel("功率密度(dB)");
ylabel("多普勒频移{\phi}(Hz)");

subplot(224)
p_map = p/max(p);
p_log = 10*log10(p_map);
plot(p_log);
title("$10log_{10}(\int s(\tau,\varphi)d\varphi)$",'interpreter','latex');
xlabel("功率密度(dB)");
ylabel("多普勒频移{\phi}(Hz)");

%% calculate CW signal in time-varible channel
dt = 128/1000/2048;
fs = 1/dt;
signal_width = 0.05;
signal_t = 0:dt:signal_width;
f_signal_CW = 1000;
rec_t = 0:dt:(signal_width+delay_range(end)/1000);
corr_t = 0:dt:(2*signal_width+delay_range(end)/1000);

CW_signal = exp(-1j*2*pi*f_signal_CW*signal_t);
rec_signal_CW = zeros(size_h(1),length(signal_t)+size_h(2)-1);
for i = 1:size_h(1)
   rec_signal_CW(i,:) = abs(conv(h(i,:),CW_signal));
end
rec_signal_CW = 20*log10(rec_signal_CW/max(rec_signal_CW,[],'all'));

for i = 1:size_h(1)
   corr_CW = xcorr(rec_signal_CW(i,:),CW_signal);
end

%% Draw signal in time-varible channel
subplot(221)
plot(signal_t,CW_signal);
subplot(222)
surf(rec_t,t_range,rec_signal_CW);
shading interp
subplot(212)
surf(corr_t,t_range,corr_CW);