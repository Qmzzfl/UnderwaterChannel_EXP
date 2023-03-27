clear;clc;close all;

%% %%%%%%%%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%
% sonar location setting
d_sonar_tang = 2;
d_sonar_normal = 2; 
d_sonar_radial = 2;

% Infinity space configure
d_infinity = 1500;

% multi-channel simulation setting
f_end_transfer = 50000;                                     %the end value of f
df_transfer = 0.01;                                         %delta f
f_transfer = 0:df_transfer:f_end_transfer-df_transfer;      %argument f of transfer function
T_transfer = f_transfer/df_transfer/f_end_transfer;         %argument t of system function(calculate according to DSP theory) 
global dt;
dt = 1/f_end_transfer;

%  multi-channel configure
H_multichannel = 100;               %deep of shallow sea
d1_multichannel = 55;               %distance between surface and source
d2_multichannel = 60;               %distance between surface and hydrophone
d_x_multichannel = 150;             %horizontal distance
zw = 1.5e6;                         %resistance of water
zb = 3.6e6;                         %resistance of seabed
global c_w;                          
c_w = 1.5e3;                       %speed of sound in water
c_b = 3e3;                          %speed of sound in seabed

%% %%%%%%%%%%%%% Signal Setting %%%%%%%%%%%%%%%%%%%%
t_signal_end = 0.1;
t_signal = 0:dt:t_signal_end;
% CW signal
f_CW = 1000;
CW = cos(2*pi*f_CW*t_signal);
CW_flip = flip(CW);
% LFM signal
f0_LFM = 1000;
f1_LFM = 1100;
LFM = chirp(t_signal,f0_LFM,t_signal_end,f1_LFM);
LFM_flip = flip(LFM);

%%%%%%%%%%%%%%%%%%%% INFINITY SPACE %%%%%%%%%%%%%%%%%%%%
% noise setting
SNR = 20;
% time list
t_observe_end = 5;
t_observe = 0:dt:t_observe_end;
t_corr = -t_observe_end:dt:t_observe_end;
% recieve list
y_rec_CW_infty1 = zeros(1,length(t_observe));
y_rec_CW_infty2 = zeros(1,length(t_observe));
y_rec_LFM_infty1 = zeros(1,length(t_observe));
y_rec_LFM_infty2 = zeros(1,length(t_observe));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%   Tangental     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_infty_tang_x = sqrt(d_infinity^2+d_sonar_tang^2);
%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve CW in infinity space when sonar place in tangental direction
[y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp] = ...
receive_infty(y_rec_CW_infty1,y_rec_CW_infty2,d_infty_tang_x,d_infty_tang_x,CW);
% Draw
figure(1)
suptitle("无限大空间中被动声呐水听器切向排布时接收信号与其互相关(CW脉冲)");
draw_signal(y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp,t_observe,t_corr);
%%%%%%%%%%%%%%%%%%%%%%%%%% LFM %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in infinity space when sonar place in tangental direction
[y_rec_LFM_infty1,y_rec_LFM_infty2,y_LFM_infty_envp] = ...
receive_infty(y_rec_LFM_infty1,y_rec_LFM_infty2,d_infty_tang_x,d_infty_tang_x,LFM);
% Draw
figure(2)
suptitle("无限大空间中被动声呐水听器切向排布时接收信号与其互相关(LFM脉冲)");
draw_signal(y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp,t_observe,t_corr);


function [signal1,signal2,signal_xcorr] = receive_infty(empty1,empty2,d1,d2,signal)
    global c_w;
    global dt;
    n_transmission1 = round(d1/c_w/dt);
    empty1(n_transmission1:n_transmission1+length(signal)-1) = signal/d1;
    signal1 = awgn(empty1,SNR,'measured');
    n_transmission2 = round(d2/c_w/dt);
    empty2(n_transmission2:n_transmission2+length(signal)-1) = signal/d2;
    signal2 = awgn(empty2,SNR,'measured');
    y_conv = xcorr(signal1,signal2);
    signal_xcorr = envelope(y_conv,30,'rms');
end

function [] = draw_signal(signal1,signal2,signal_xcorr,t_signal,t_xcorr)
    subplot('position',[0.13,0.67,0.8,0.2])
    plot(t_signal,signal1);
    xlabel("时间/s");
    ylabel("信号幅值");
    title("水听器1接收信号幅值");
    subplot('position',[0.13,0.38,0.8,0.2])
    plot(t_signal,signal2);
    xlabel("时间/s");
    ylabel("信号幅值");
    title("水听器2接收信号幅值");
    subplot('position',[0.13,0.09,0.8,0.2])
    plot(t_xcorr,mapminmax(signal_xcorr,0));
    xlim([t_xcorr(1) t_xcorr(end)]);
    xlabel("两接收波形时延/s");
    ylabel("互相关归一化幅值");
    title("两接收信号互相关结果");
end