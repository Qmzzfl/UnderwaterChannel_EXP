clear;clc;close all;

H = 100;                %deep of shallow sea
d1 = 55;                %distance between surface and source
d2 = 60;                 %distance between surface and hydrophone
d_x = 1000;                %horizontal distance
zw = 1.5e6;             %resistance of water
zb = 3.6e6;              %resistance of seabed
%Rp = (zb-zw)/(zb+zw);   %reflect coefficient(wrong!oblique incidence!)
cw = 1.5e3;              %velocity of sound
cb = 3e3;

f_end = 30000;           %the end value of f
df = 0.01;              %delta f
f = 0:df:f_end;         %argument f of transfer function
T = f/df/f_end;         %argument t of system function(calculate according to DSP theory) 




%% %%%%%%%%%%%%% SINGLE EXP %%%%%%%%%%%%%%%%%%%%%%
f_observeRange = [0,1000];
t_observeRange = [0.6,2];
ref_num = 10;
Hf = Transfer(ref_num,f,d1,d2,d_x,H,cw,cb,zw,zb);
ht = ifft(Hf,'symmetric'); 
%% Draw single exp
figure(1)
subplot('position',[0.1 0.757 0.8 0.2])
% plot(f,20*log(abs(Hf)));
plot(f,abs(Hf));
xlim(f_observeRange)
title("多途信道的幅频响应");
xlabel("Frequency/Hz");
ylabel("Relative Magnitude/dB");
subplot('position',[0.1 0.427 0.8 0.2])
Hf_angle = angle(Hf)/pi*180;
plot(f,Hf_angle);
xlim(f_observeRange)
yticks([-180 -90 0 90 180])
title("多途信道的相频特性");
xlabel("Frequency/Hz");
ylabel("Phase/degree");
subplot('position',[0.1 0.07 0.8 0.2])
plot(T,ht);
xlim(t_observeRange)
title("多途信道的冲激响应函数");
xlabel("Time/s");
ylabel("Relative Magnitude");

%%  %%%%%%%%%%%% COMPARE EXP %%%%%%%%%%%%%%%%%%%%%
f_observeRange = [0,1000];
ref_num = 5;
%% %%%%%%%%%% the depth of SOURCE %%%%%%%%%%%%%%%
figure(2)
compare_range = 0:0.05:0.5;
hf_compare = zeros(round(length(compare_range)),length(f));
j = 1;
for i = compare_range
    hf_compare(j,:) = abs(Transfer(ref_num,f,d1+i,d2,d_x,H,cw,cb,zw,zb));
    j = j + 1;
end

for i = 1:length(compare_range)
    subplot('position',[0.1 0.97-i*0.076 0.8 0.075])
    plot(f,hf_compare(i,:));
    xlim(f_observeRange)
    if i ~= length(compare_range)
        set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    else 
        set(gca,'ytick',[],'yticklabel',[])
    end
end
suptitle('传输函数随换能器深度增加的变化情况');
xlabel("Frequency/Hz");

figure(10)
subplot(411)
compare_similar = zeros(1,length(compare_range));
for i = 1:length(compare_range)
    compare_matrix = corrcoef(hf_compare(1,:),hf_compare(i,:));
    compare_similar(1,i) = abs(compare_matrix(1,2));
end
plot(compare_range(2:end),compare_similar(2:end));
title("改变发射器深度后的传输函数相似程度");
xlabel("变化量/m");
ylabel("相似度");
ylim([0 0.6]);


%% %%%%%%%%%% the depth of HYDROPHONE %%%%%%%%%%%%%%%
figure(3)
compare_range = 0:0.05:0.5;
hf_compare = zeros(round(length(compare_range)),length(f));
j = 1;
for i = compare_range
    hf_compare(j,:) = abs(Transfer(ref_num,f,d1,d2+i,d_x,H,cw,cb,zw,zb));
    j = j + 1;
end

for i = 1:length(compare_range)
    subplot('position',[0.1 0.97-i*0.075 0.8 0.075])
    plot(f,hf_compare(i,:));
    xlim(f_observeRange)
    if i ~= length(compare_range)
        set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    else 
        set(gca,'ytick',[],'yticklabel',[])
    end
end
suptitle('传输函数随接收器深度增加的变化情况');
xlabel("Frequency/Hz");

figure(10)
subplot(412)
compare_similar = zeros(1,length(compare_range));
for i = 1:length(compare_range)
    compare_matrix = corrcoef(hf_compare(1,:),hf_compare(i,:));
    compare_similar(1,i) = abs(compare_matrix(1,2));
end
plot(compare_range(2:end),compare_similar(2:end));
title("改变接收器深度后的传输函数相似程度");
xlabel("变化量/m");
ylabel("相似度");
ylim([0 0.6]);

%% %%%%%%%%%% the depth of SEA %%%%%%%%%%%%%%%
figure(4)
compare_range = 0:0.05:0.5;
hf_compare = zeros(round(length(compare_range)),length(f));
j = 1;
for i = compare_range
    hf_compare(j,:) = abs(Transfer(ref_num,f,d1,d2,d_x,H+i,cw,cb,zw,zb));
    j = j + 1;
end

for i = 1:length(compare_range)
    subplot('position',[0.1 0.97-i*0.075 0.8 0.075])
    plot(f,hf_compare(i,:));
    xlim(f_observeRange)
    if i ~= length(compare_range)
        set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    else 
        set(gca,'ytick',[],'yticklabel',[])
    end
end
suptitle('传输函数随海洋深度增加的变化情况');
xlabel("Frequency/Hz");

figure(10)
subplot(413)
compare_similar = zeros(1,length(compare_range));
for i = 1:length(compare_range)
    compare_matrix = corrcoef(hf_compare(1,:),hf_compare(i,:));
    compare_similar(1,i) = abs(compare_matrix(1,2));
end
plot(compare_range(2:end),compare_similar(2:end));
title("改变海洋深度后的传输函数相似程度");
xlabel("变化量/m");
ylabel("相似度");
ylim([0 0.6]);

%% %%%%%%%%%% the Horizontal DISTANCE %%%%%%%%%%%%%%%
figure(5)
compare_range = 0:1:10;
hf_compare = zeros(round(length(compare_range)),length(f));
j = 1;
for i = compare_range
    hf_compare(j,:) = abs(Transfer(ref_num,f,d1,d2,d_x+i,H,cw,cb,zw,zb));
    j = j + 1;
end

for i = 1:length(compare_range)
    subplot('position',[0.1 0.97-i*0.075 0.8 0.075])
    plot(f,hf_compare(i,:));
    xlim(f_observeRange)
    if i ~= length(compare_range)
        set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    else 
        set(gca,'ytick',[],'yticklabel',[])
    end
end
suptitle('传输函数随水平距离增加的变化情况');
xlabel("Frequency/Hz");

figure(10)
subplot(414)
compare_similar = zeros(1,length(compare_range));
for i = 1:length(compare_range)
    compare_matrix = corrcoef(hf_compare(1,:),hf_compare(i,:));
    compare_similar(1,i) = abs(compare_matrix(1,2));
end
plot(compare_range(2:end),compare_similar(2:end));
title("改变水平距离后的传输函数相似程度");
xlabel("变化量/m");
ylabel("相似度");
ylim([0 0.6]);


% %%  %%%%%%%%%%%% COMPARE EXP %%%%%%%%%%%%%%%%%%%%%
% compare_range = 0:1:5;
% ht = zeros(round(length(compare_range)),length(T));
% ref_num = 5;
% %%%%%%%%%%%% the depth of SOURCE %%%%%%%%%%%%%%%
% figure(2)
% j = 1;
% for i = compare_range
%     Hf = Transfer(ref_num,f,d1+i,d2,d_x,H,cw,cb,zw,zb);
%     ht(j,:) = ifft(Hf,'symmetric'); 
%     j = j + 1;
% end
% 
% for i = 1:length(compare_range)
%     count = 1;
%     u = zeros(length(T),3);
%     for j = T
%         u(count,1) = i;
%         u(count,2) = j;
%         u(count,3) = ht(i,count);
%         count = count + 1;
%     end    
%     x = u(:,1);
%     y = u(:,2);
%     z = u(:,3);
%     plot3(x,y,z);
%     ylim(t_observeRange);
%     grid on;
%     hold on;
% end
% view(90,80)
% 
% %%%%%%%%%%%% the depth of HYDROPHONE %%%%%%%%%%%%%%%
% figure(3)
% j = 1;
% for i = compare_range
%     Hf = Transfer(ref_num,f,d1,d2+i,d_x,H,cw,cb,zw,zb);
%     ht(j,:) = ifft(Hf,'symmetric'); 
%     j = j + 1;
% end
% 
% for i = 1:length(compare_range)
%     count = 1;
%     u = zeros(length(T),3);
%     for j = T
%         u(count,1) = i;
%         u(count,2) = j;
%         u(count,3) = ht(i,count);
%         count = count + 1;
%     end    
%     x = u(:,1);
%     y = u(:,2);
%     z = u(:,3);
%     plot3(x,y,z);
%     ylim(t_observeRange);
%     grid on;
%     hold on;
% end
% view(90,80)
% 
% %%%%%%%%%%%% the depth of SEA %%%%%%%%%%%%%%%
% figure(4)
% j = 1;
% for i = compare_range
%     Hf = Transfer(ref_num,f,d1,d2,d_x,H+i,cw,cb,zw,zb);
%     ht(j,:) = ifft(Hf,'symmetric'); 
%     j = j + 1;
% end
% 
% for i = 1:length(compare_range)
%     count = 1;
%     u = zeros(length(T),3);
%     for j = T
%         u(count,1) = i;
%         u(count,2) = j;
%         u(count,3) = ht(i,count);
%         count = count + 1;
%     end    
%     x = u(:,1);
%     y = u(:,2);
%     z = u(:,3);
%     plot3(x,y,z);
%     ylim(t_observeRange);
%     grid on;
%     hold on;
% end
% view(90,80)
% 
% %%%%%%%%%%%% the Horizontal DISTANCE %%%%%%%%%%%%%%%
% figure(5)
% j = 1;
% for i = compare_range
%     Hf = Transfer(ref_num,f,d1,d2,d_x+i,H,cw,cb,zw,zb);
%     ht(j,:) = ifft(Hf,'symmetric'); 
%     j = j + 1;
% end
% 
% for i = 1:length(compare_range)
%     count = 1;
%     u = zeros(length(T),3);
%     for j = T
%         u(count,1) = i;
%         u(count,2) = j;
%         u(count,3) = ht(i,count);
%         count = count + 1;
%     end    
%     x = u(:,1);
%     y = u(:,2);
%     z = u(:,3);
%     plot3(x,y,z);
%     ylim(t_observeRange);
%     grid on;
%     hold on;
% end
% view(90,80)




%%%%%%%%%% Function Definition %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name: SurfaceReflect
%Function parameter:the coordinate of source before reflect
%Function return:the coordinate of vitaul source after reflect by SURFACE
function z_SurfaceReflect = SurfaceReflect(z_origin)
    z_SurfaceReflect = -z_origin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name: SurfaceReflect
%Function parameter:the coordinate of source before reflect
%Function return:the coordinate of vitaul source after reflect by SEABED
function z_SeabedReflect = SeabedReflect(z_origin,Depth)
    z_SeabedReflect = 2*Depth-z_origin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function name:Transfer
%Output the transfer function of the multipath channel correspond by the
%parameters
%Function parameter:
%   reflection_times:reflection times of channel
%   f_list:the 1-D matrix of argument f
%   z_source:the coordinate of source
%   z_hydrophone:the coordinate of source
%   horizon_d: horizon distance between source and hydrophone
%   Depth:the depth of sea
%   c:wave speed in media 
%   z1-media resistance  z2-seabed resistance
%   c1-media sound speed c2-seabed sound speed
%Function Return:The transfer function of the multichannel
function Transfer_Function = Transfer(reflection_times,...
                                      f_list,...
                                      z_source,z_hydrophone,...
                                      horizon_d,...
                                      Depth,c1,c2,z1,z2)
    %calculate the distance of direcet wave
    Distance_direct = sqrt((z_source-z_hydrophone)^2+horizon_d^2);
    %Define matrix of coordinate of vitual source and correspond magnitude
    z_virtualSource = zeros(2,reflection_times);
    mag_virtualSource = zeros(2,reflection_times);
    Distance = zeros(2,reflection_times);
    %Define the one-reflection parameters
    z_virtualSource(1,1) = SurfaceReflect(z_source);
    mag_virtualSource(1,1) = -1; 
    z_virtualSource(2,1) = SeabedReflect(z_source,Depth);
    mag_virtualSource(2,1) = 1; 
    % cross operation to calculate other parameters 
    for i = 2:reflection_times
        z_virtualSource(1,i) = SurfaceReflect(z_virtualSource(2,i-1));
        Distance(1,i) = sqrt((z_hydrophone-z_virtualSource(1,i)).^2+horizon_d^2);
        mag_virtualSource(1,i) = -1*mag_virtualSource(2,i-1); 

        z_virtualSource(2,i) = SeabedReflect(z_virtualSource(1,i-1),Depth);
        Distance(2,i) = sqrt((z_hydrophone-z_virtualSource(2,i)).^2+horizon_d^2);
        mag_virtualSource(2,i) = (refCo(Distance(2,i),horizon_d,z1,z2,c1,c2))^floor(i/2);
        %mag_virtualSource(2,i) = -1*mag_virtualSource(1,i-1); 
    end
    %Calculate the Distance between vitual sources and hyrdophone
    Distance = sqrt((z_hydrophone-z_virtualSource).^2+horizon_d^2);
    %Put it into 1-D for easy indexing
    Distance = reshape(Distance,1,[]);
    mag_virtualSource = reshape(mag_virtualSource,1,[]);
    
    %Accumulate to get transfer function
    Transfer_Function = 1/Distance_direct*...
            exp(-1j*2*pi*f_list*Distance_direct/c1);
    for i = 1:length(Distance)
        Transfer_Function = Transfer_Function + ...
            mag_virtualSource(i)/Distance(i)*...
            exp(-1j*2*pi*f_list*Distance(i)/c1);
    end
end

% Function name:refCo
% Work:Calculate Reflect Coefficient in different channel
% Input:L-channel length  x-Horizontal distance  
%       z1-media resistance  z2-seabed resistance
%       c1-media sound speed c2-seabed sound speed
% Return:reflect coefficient
function Rp = refCo(L,x,z1,z2,c1,c2)
    theta_i = asin(x/L);
    theta_t = asin(c2*x/(c1*L));
    z1_ob = z1/cos(theta_i);
    z2_ob = z2/cos(theta_t);
    Rp = real((z2_ob-z1_ob)/(z2_ob+z1_ob));
end

