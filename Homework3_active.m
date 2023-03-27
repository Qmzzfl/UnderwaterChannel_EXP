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

%% %%%%%%%%%%%%% INFINITY-CHANNEL %%%%%%%%%%%%%%%%%%%
SNR = -30;
% time list
t_observe_end = 5;
t_observe = 0:dt:t_observe_end;
t_corr = 0:dt:t_observe_end+t_signal_end;
% recieve list
y_rec_CW_infty = zeros(1,length(t_observe));
y_rec_LFM_infty = zeros(1,length(t_observe));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve CW in infinity space
n_transmission_CW_infty = round(d_infinity/c_w/dt);
y_rec_CW_infty(n_transmission_CW_infty:n_transmission_CW_infty+length(CW)-1) = CW/d_infinity;
y_rec_CW_infty = awgn(y_rec_CW_infty,SNR,'measured');
y_corr_CW_infty = conv(CW_flip,y_rec_CW_infty);
% Draw
figure(1)
suptitle("CW���������޴�ռ��еĴ���");
subplot('position',[0.13,0.67,0.8,0.2])
plot(t_signal,CW);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.38,0.8,0.2])
plot(t_observe,y_rec_CW_infty);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.09,0.8,0.2])
plot(t_corr,y_corr_CW_infty);
xlim([0 t_observe_end+t_signal_end]);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źž��������������Ĳ���");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in infinity space
n_transmission_LFM_infty = round(d_infinity/c_w/dt);
y_rec_LFM_infty(n_transmission_LFM_infty:n_transmission_LFM_infty+length(LFM)-1) = LFM/d_infinity;
y_rec_LFM_infty = awgn(y_rec_LFM_infty,SNR,'measured');
y_corr_LFM_infty = conv(LFM_flip,y_rec_LFM_infty);
% Draw
figure(2)
suptitle("LFM���������޴�ռ��еĴ���");
subplot('position',[0.13,0.67,0.8,0.2])
plot(t_signal,LFM);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.38,0.8,0.2])
plot(t_observe,y_rec_LFM_infty);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.09,0.8,0.2])
plot(t_corr,y_corr_LFM_infty);
xlim([0 t_observe_end+t_signal_end]);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źž��������������Ĳ���");

%% %%%%%%%%%%%%% MULTI-CHANNEL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-channel setting and calculation
SNR = -30;
t_observe_end = 5;
t_observe = 0:dt:t_observe_end;
t_ht_end = 10;
n_ht_end = t_ht_end*5000000/100;
t_conv = 0:dt:t_ht_end+t_signal_end-dt;
t_corr = 0:dt:t_conv(end)++t_signal_end;
ref_num = 3;
Hf = Transfer(ref_num,f_transfer,d1_multichannel,d2_multichannel,d_x_multichannel,...
                H_multichannel,c_w,c_b,zw,zb);
ht = ifft(Hf,'symmetric'); 
% Draw multi-channel
figure(3)
suptitle("��;�ŵ�������");
subplot('position',[0.13,0.51,0.8,0.32])
plot(f_transfer,abs(Hf));
xlim([0 f_end_transfer]);
xlabel("Ƶ��/Hz");
ylabel("���亯����ֵ");
title("��;�ŵ��ķ�Ƶ��Ӧ");
subplot('position',[0.13,0.08,0.8,0.32])
plot(T_transfer,ht);
xlim([0 t_observe_end]);
xlabel("ʱ��/s");
ylabel("�弤��Ӧ������ֵ");
title("��;�ŵ��ĳ弤��Ӧ����");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve CW in multi-channel 
y_rec_CW_multichannel = conv(ht(1:n_ht_end),CW);
y_rec_CW_multichannel = awgn(y_rec_CW_multichannel,SNR,'measured');
y_corr_CW_multichannel = conv(CW_flip,y_rec_CW_multichannel);
% Draw
figure(4)
suptitle("CW��������ɶ�;�ŵ��еĴ���");
subplot('position',[0.13,0.67,0.8,0.2])
plot(t_signal,CW);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.38,0.8,0.2])
plot(t_conv,y_rec_CW_multichannel);
xlim([0 t_conv(end)]);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.09,0.8,0.2])
plot(t_corr,y_corr_CW_multichannel);
xlim([0 t_observe_end+t_signal_end]);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źž��������������Ĳ���");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve LFM in multi-channel
y_rec_LFM_multichannel = conv(ht(1:500000),LFM);
y_rec_LFM_multichannel = awgn(y_rec_LFM_multichannel,SNR,'measured');
y_corr_LFM_multichannel = conv(LFM_flip,y_rec_LFM_multichannel);
% Draw
figure(5)
suptitle("LFM��������ɶ�;�ŵ��еĴ���");
subplot('position',[0.13,0.67,0.8,0.2])
plot(t_signal,LFM);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.38,0.8,0.2])
plot(t_conv,y_rec_LFM_multichannel);
xlim([0 t_conv(end)]);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źŲ���");
subplot('position',[0.13,0.09,0.8,0.2])
plot(t_corr,y_corr_LFM_multichannel);
xlim([0 t_observe_end+t_signal_end]);
xlabel("ʱ��/s");
ylabel("�źŷ�ֵ");
title("�����źž��������������Ĳ���");