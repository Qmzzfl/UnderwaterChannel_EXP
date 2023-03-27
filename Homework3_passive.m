clear;clc;close all;

% Infinity space configure
d_infinity = 1500;

%  multi-channel configure
H_multichannel = 100;               %deep of shallow sea
d1_multichannel = 55;               %distance between surface and source
d2_multichannel = 60;               %distance between surface and hydrophone
d_x_multichannel = 150;             %horizontal distance
zw = 1.5e6;                         %resistance of water
zb = 3.6e6;                         %resistance of seabed
c_w = 1.5e3;                         %speed of sound in water
c_b = 3e3;                           %speed of sound in seabed

% multi-channel simulation setting
f_end_transfer = 50000;                                     %the end value of f
df_transfer = 0.01;                                         %delta f
f_transfer = 0:df_transfer:f_end_transfer-df_transfer;      %argument f of transfer function
T_transfer = f_transfer/df_transfer/f_end_transfer;         %argument t of system function(calculate according to DSP theory) 
dt = 1/f_end_transfer;

% sonar location setting
d_sonar_tang = 2;
d_sonar_normal = 2; 
d_sonar_radial = 2;

%% %%%%%%%%%%%%% Signal Setting %%%%%%%%%%%%%%%%%%%%
t_signal_end = 0.1;
t_signal = 0:dt:t_signal_end;
% CW signal
f_CW = 1000;
CW = cos(2*pi*f_CW*t_signal);
CW_flip = flip(CW);
% LFM signal
f0_LFM = 1000;
f1_LFM = 2000;
LFM = chirp(t_signal,f0_LFM,t_signal_end,f1_LFM);
LFM_flip = flip(LFM);

%%%%%%%%%%%%%%%%%%%% INFINITY SPACE %%%%%%%%%%%%%%%%%%%%
% noise setting
SNR = -20;
% time list
t_observe_end = 15;
t_observe = 0:dt:t_observe_end;
t_corr = 0:dt:2*t_observe_end;
% recieve list
y_rec_CW_infty1 = zeros(1,length(t_observe));
y_rec_CW_infty2 = zeros(1,length(t_observe));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%   Tangental     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_infty_tang_x = sqrt(d_infinity^2+d_sonar_tang^2);
%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve CW in infinity space when sonar place in tangental place
n_transmission_CW_infty = round(d_infty_tang_x/c_w/dt);
y_rec_CW_infty1(n_transmission_CW_infty:n_transmission_CW_infty+length(CW)-1) = CW/d_infinity;
y_rec_CW_infty1 = awgn(y_rec_CW_infty1,SNR,'measured');
y_rec_CW_infty2(n_transmission_CW_infty:n_transmission_CW_infty+length(CW)-1) = CW/d_infinity;
y_rec_CW_infty2 = awgn(y_rec_CW_infty2,SNR,'measured');
y_CW_infty_conv = conv(y_rec_CW_infty1,flip(y_rec_CW_infty2));
% Draw
figure(1)
suptitle("无限大空间中被动声呐水听器切向排布时接收信号与其互相关");
subplot('position',[0.13,0.67,0.8,0.2])
plot(t_observe,y_rec_CW_infty1);
xlabel("时间/s");
ylabel("信号幅值");
title("发射信号波形");
subplot('position',[0.13,0.38,0.8,0.2])
plot(t_observe,y_rec_CW_infty2);
xlabel("时间/s");
ylabel("信号幅值");
title("接收信号波形");
subplot('position',[0.13,0.09,0.8,0.2])
plot(t_corr,y_CW_infty_conv);
xlim([0 t_corr(end)]);
xlabel("时间/s");
ylabel("信号幅值");
title("接收信号经过拷贝相关器后的波形");