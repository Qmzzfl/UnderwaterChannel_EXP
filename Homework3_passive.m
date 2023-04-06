clear;clc;close all;

%% %%%%%%%%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%
% sonar location setting
d_sonar_tang = 2;
d_sonar_normal = 150; 
d_sonar_radial = 150;

% Infinity space configure
d_infinity_y = 500;
d_infinity_z = 100;

% multi-channel simulation setting
f_end_transfer = 5000;                                     %the end value of f
df_transfer = 0.1;                                         %delta f
f_transfer = 0:df_transfer:f_end_transfer-df_transfer;      %argument f of transfer function
T_transfer = f_transfer/df_transfer/f_end_transfer;         %argument t of system function(calculate according to DSP theory) 
global dt;
dt = 1/f_end_transfer;

%  multi-channel configure
H_multichannel = 100;               %deep of shallow sea
d1_multichannel = 50;               %distance between surface and source
d2_multichannel = 50;               %distance between surface and hydrophone
d_x_multichannel = 500;             %horizontal distance
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
% noise setting
global SNR;
SNR = 40;

%% %%%%%%%%%%%%%%%%%% INFINITY SPACE %%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%   Tangental     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_infty_tang_x = sqrt(d_infinity_y^2+d_infinity_z^2+d_sonar_tang^2);
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
draw_signal(y_rec_LFM_infty1,y_rec_LFM_infty2,y_LFM_infty_envp,t_observe,t_corr);

%%%%%%%%%%%%%%%%   Normal         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_infty_norm_x1 = sqrt(d_infinity_y^2+(d_infinity_z-d_sonar_normal)^2);
d_infty_norm_x2 = sqrt(d_infinity_y^2+(d_infinity_z+d_sonar_normal)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve CW in infinity space when sonar place in tangental direction
[y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp] = ...
receive_infty(y_rec_CW_infty1,y_rec_CW_infty2,d_infty_norm_x1,d_infty_norm_x2,CW);
% Draw
figure(3)
suptitle("无限大空间中被动声呐水听器法向排布时接收信号与其互相关(CW脉冲)");
draw_signal(y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp,t_observe,t_corr);
%%%%%%%%%%%%%%%%%%%%%%%%%% LFM %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in infinity space when sonar place in tangental direction
[y_rec_LFM_infty1,y_rec_LFM_infty2,y_LFM_infty_envp] = ...
receive_infty(y_rec_LFM_infty1,y_rec_LFM_infty2,d_infty_norm_x1,d_infty_norm_x2,LFM);
% Draw
figure(4)
suptitle("无限大空间中被动声呐水听器法向排布时接收信号与其互相关(LFM脉冲)");
draw_signal(y_rec_LFM_infty1,y_rec_LFM_infty2,y_LFM_infty_envp,t_observe,t_corr);

%%%%%%%%%%%%%%%%   radioal         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_infty_rad_x1 = sqrt(d_infinity_z^2+(d_infinity_y-d_sonar_radial)^2);
d_infty_rad_x2 = sqrt(d_infinity_z^2+(d_infinity_y+d_sonar_radial)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve CW in infinity space when sonar place in tangental direction
[y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp] = ...
receive_infty(y_rec_CW_infty1,y_rec_CW_infty2,d_infty_rad_x1,d_infty_rad_x2,CW);
% Draw
figure(5)
suptitle("无限大空间中被动声呐水听器径向排布时接收信号与其互相关(CW脉冲)");
draw_signal(y_rec_CW_infty1,y_rec_CW_infty2,y_CW_infty_envp,t_observe,t_corr);
%%%%%%%%%%%%%%%%%%%%%%%%%% LFM %%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in infinity space when sonar place in tangental direction
[y_rec_LFM_infty1,y_rec_LFM_infty2,y_LFM_infty_envp] = ...
receive_infty(y_rec_LFM_infty1,y_rec_LFM_infty2,d_infty_rad_x1,d_infty_rad_x2,LFM);
% Draw
figure(6)
suptitle("无限大空间中被动声呐水听器径向排布时接收信号与其互相关(LFM脉冲)");
draw_signal(y_rec_LFM_infty1,y_rec_LFM_infty2,y_LFM_infty_envp,t_observe,t_corr);

%% %%%%%%%%%%%%%%%%%%%%%%% MULTI-CHANNEL %%%%%%%%%%%%%%%%%%%%%
% Multi-channel setting
t_observe_end = 5;
t_observe = 0:dt:t_observe_end;
t_ht_end = 5;
n_ht_end = t_ht_end*length(T_transfer)*(df_transfer);
t_conv = 0:dt:t_ht_end+t_signal_end-dt;
t_corr = -t_conv(end):dt:t_conv(end);
ref_num = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Tangental     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_x_multichannel_tang = sqrt(d_x_multichannel^2+d_sonar_tang^2);
ht1 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel_tang,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
ht2 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel_tang,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
%%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multichannel when sonar place in tangental direction
[y_rec_CW_multichannel1,y_rec_CW_multichannel2,y_CW_multichannel_envp] = ...
receive_multichannel(ht1(1:n_ht_end),ht2(1:n_ht_end),CW);
% Draw
figure(7)
suptitle("相干多途信道中被动声呐水听器切向排布时接收信号与其互相关(CW脉冲)")
draw_signal(y_rec_CW_multichannel1,y_rec_CW_multichannel2,y_CW_multichannel_envp,...
            t_conv,t_corr);
%%%%%%%%%%%%%%%%%%%%%%%%%%% LFM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multichannel when sonar place in tangental direction
[y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp] = ...
receive_multichannel(ht1(1:n_ht_end),ht2(1:n_ht_end),LFM);
% Draw
figure(8)
suptitle("相干多途信道中被动声呐水听器切向排布时接收信号与其互相关(LFM脉冲)")
draw_signal(y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp,...
            t_conv,t_corr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   Normal     %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ht1 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel-d_sonar_normal,...
                    d_x_multichannel,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
ht2 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel+d_sonar_normal,...
                    d_x_multichannel,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
%%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multichannel when sonar place in tangental direction
[y_rec_CW_multichannel1,y_rec_CW_multichannel2,y_CW_multichannel_envp] = ...
receive_multichannel(ht1(1:n_ht_end),ht2(1:n_ht_end),CW); 
% Draw
figure(9)
suptitle("相干多途信道中被动声呐水听器法向排布时接收信号与其互相关(CW脉冲)")
draw_signal(y_rec_CW_multichannel1,y_rec_CW_multichannel2,y_CW_multichannel_envp,...
            t_conv,t_corr);
%%%%%%%%%%%%%%%%%%%%%%%%%%% LFM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multichannel when sonar place in tangental direction
[y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp] = ...
receive_multichannel(ht1(1:n_ht_end),ht2(1:n_ht_end),LFM);
% Draw
figure(10)
suptitle("相干多途信道中被动声呐水听器法向排布时接收信号与其互相关(LFM脉冲)")
draw_signal(y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp,...
            t_conv,t_corr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   radial         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ht1 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel-d_sonar_radial,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
ht2 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel+d_sonar_radial,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
%%%%%%%%%%%%%%%%%%%%%%%%%%% CW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multichannel when sonar place in tangental direction
[y_rec_CW_multichannel1,y_rec_CW_multichannel2,y_CW_multichannel_envp] = ...
receive_multichannel(ht1(1:n_ht_end),ht2(1:n_ht_end),CW);
% Draw
figure(11)
suptitle("相干多途信道中被动声呐水听器径向排布时接收信号与其互相关(CW脉冲)")
draw_signal(y_rec_CW_multichannel1,y_rec_CW_multichannel2,y_CW_multichannel_envp,...
            t_conv,t_corr);
%%%%%%%%%%%%%%%%%%%%%%%%%%% LFM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multichannel when sonar place in tangental direction
[y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp] = ...
receive_multichannel(ht1(1:n_ht_end),ht2(1:n_ht_end),LFM);
% Draw
figure(12)
suptitle("相干多途信道中被动声呐水听器径向排布时接收信号与其互相关(LFM脉冲)")
draw_signal(y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp,...
            t_conv,t_corr);

%% multi-channel multi exp
figure(13)
test_domain = 0:0.1:20;
j = 1;
test_domain_len = length(test_domain);
test_res = zeros(1,test_domain_len);

% radial
for i = test_domain
    ht_test1 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel-i,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
    ht_test2 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel+i,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
   [y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp] = ...
        receive_multichannel(ht_test1(1:n_ht_end),ht_test2(1:n_ht_end),LFM);
%     corr_test = corrcoef(y_rec_LFM_multichannel1,y_rec_LFM_multichannel2);
%     test_res(j) = abs(corr_test(1,2));
    test_res(j) = max(y_LFM_multichannel_envp);
    j = j + 1;
end
hold on
plot(test_domain,test_res);

% normal
ht_test1 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                d_x_multichannel,H_multichannel,c_w,c_b,zw,zb),...
       'symmetric');
j = 1;
for i = test_domain

    ht_test2 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel+i,...
                    d_x_multichannel,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
   [y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp] = ...
        receive_multichannel(ht_test1(1:n_ht_end),ht_test2(1:n_ht_end),LFM);
%     corr_test = corrcoef(y_rec_LFM_multichannel1,y_rec_LFM_multichannel2);
%     test_res(j) = abs(corr_test(1,2));
    test_res(j) = max(y_LFM_multichannel_envp);
    j = j + 1;
end
plot(test_domain,test_res);

% tangental
ht_test1 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                d_x_multichannel,H_multichannel,c_w,c_b,zw,zb),...
       'symmetric');
j = 1;
for i = test_domain
    d_x_multichannel_tang = sqrt(d_x_multichannel^2+i^2);
    ht_test2 = ifft(Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,...
                    d_x_multichannel_tang,H_multichannel,c_w,c_b,zw,zb),...
           'symmetric');
   [y_rec_LFM_multichannel1,y_rec_LFM_multichannel2,y_LFM_multichannel_envp] = ...
        receive_multichannel(ht_test1(1:n_ht_end),ht_test2(1:n_ht_end),LFM);
%     corr_test = corrcoef(y_rec_LFM_multichannel1,y_rec_LFM_multichannel2);
%     test_res(j) = abs(corr_test(1,2));
    test_res(j) = max(y_LFM_multichannel_envp);
    j = j + 1;
end
plot(test_domain,test_res);
legend(["径向","法向","切向"])

function [signal1,signal2,signal_xcorr] = receive_infty(empty1,empty2,d1,d2,signal)
    global c_w dt SNR;
    empty1(:)=0;
    empty2(:)=0;
    n_transmission1 = round(d1/c_w/dt);
    empty1(n_transmission1:n_transmission1+length(signal)-1) = signal/d1;
    signal1 = awgn(empty1,SNR,'measured');
    n_transmission2 = round(d2/c_w/dt);
    empty2(n_transmission2:n_transmission2+length(signal)-1) = signal/d2;
    signal2 = awgn(empty2,SNR,'measured');
    y_conv = xcorr(signal1,signal2);
    signal_xcorr = envelope(y_conv,30,'rms');
end

function [signal1,signal2,signal_xcorr] = receive_multichannel(ht1,ht2,signal)
    global SNR;
    signal1 = awgn(conv(ht1,signal),SNR,'measured');
    signal2 = awgn(conv(ht2,signal),SNR,'measured');
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