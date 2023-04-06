clc;clear;close all;

t_range = linspace(0,32.9,258);
df_max = 258/32.9/2;
df_range = linspace(-df_max,df_max,258);
delay_range = linspace(0,128,2048);
fai_range = linspace(-8,8,2048);
load('.\NOF1\mat\NOF1_001.mat');

%% calculate h
figure(1)
subplot(231)
h_abs = mapminmax(abs(h),0,1);
h_log = 20*log10(h_abs);
load clown
image(delay_range,t_range,h_log(1:258,1:4:end),'CDataMapping','scaled');
colormap(slanCM(153,160))
caxis([-40 0])
% colorbar

%% calculate S
size_h = size(h);
s = zeros(size(h));
for i = 1:size_h(2)
    s(:,i) = abs(fftshift(fft(h(:,i))));
end
s_map = s/max(s,[],'all');
s_log = 20*log10(s_map);
subplot(232)
load clown
image(delay_range,df_range,s_log,'CDataMapping','scaled');
colormap(slanCM(153,160))
caxis([-40 0])
colorbar

%% calculate B
B = zeros(size_h);
for i = 1:size_h(1)
    B(i,:) = fftshift(fft(s(i,:)));
end
subplot(233)
load clown
image(fai_range,df_range,20*log10(mapminmax(abs(B),0,1)),'CDataMapping','scaled')
colormap(slanCM(153,160))
caxis([-40 0])
colorbar

Rb0_abs = abs(B(:,1025));
Rb0_map = Rb0_abs/max(Rb0_abs);
Rb0_log = 20*log10(Rb0_map);
subplot(223)
plot(Rb0_log,df_range)
%% calculate Rh
p = zeros(1,2048);
for i = 1:2048
    p(i) = (sum(s(:,i)))^2;
end
subplot(224)
p_map = p/max(p);
p_log = 10*log10(p_map);
plot(p_log);

