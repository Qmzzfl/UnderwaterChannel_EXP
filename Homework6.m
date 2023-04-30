clear;close all;clc;

% --------------------- zw --------------------
%     |   |                             |
%     |   ds                            dr
%     |   |                             |
%     H   o----------- d_x -------------o
%     |    
%     |
%-----------------------zb--------------------

H = 500;                %deep of shallow sea
ds = 250;                %distance between surface and source
dr = 250;                 %distance between surface and hydrophone
d_x = 100;                %horizontal distance
zw = 1.5e6;             %resistance of water
zb = 3.6e6;              %resistance of seabed
%Rp = (zb-zw)/(zb+zw);   %reflect coefficient(wrong!oblique incidence!)
cw = 1.5e3;              %velocity of sound
cb = 3e3;

f_end = 30000;           %the end value of f
df = 0.01;              %delta f
f = 0:df:f_end;         %argument f of transfer function
fs = 1/f_end;
T = f/df/f_end;         %argument t of system function(calculate according to DSP theory) 

signal_t_end = 0.1;
signal_t = 0:fs:signal_t_end;

CW_f = 250;
CW = sin(2*pi*CW_f*signal_t);

LFM_f_start = 100;
LFM_f_end = 1000;
LFM = chirp(signal_t,LFM_f_start,signal_t_end,LFM_f_end);

element_num = 10;
element_interval = 3;
num_center = floor(element_num/2);
element_position = linspace(-num_center*element_interval,...
                             num_center*element_interval,...
                             element_num);
%%                          
ht = zeros(1,length(T));
for i = 1:element_num
    element_dx = sqrt(d_x^2+element_position(i)^2);
    ht_temp = ifft(Transfer(10,f,ds,dr,element_dx,H,cw,cb,zw,cb),...
                'symmetric');
    ht = ht + ht_temp;
end
hold on
plot(conv(ht,CW));
%%
ht = zeros(1,length(T));
for i = 1:element_num
    ht_temp = ifft(Transfer(10,f,ds,dr+element_position(i),d_x,H,cw,cb,zw,cb),...
                'symmetric');
    ht = ht + ht_temp;
end
plot(conv(ht,CW));
%%
ht = zeros(1,length(T));
for i = 1:element_num
    ht_temp = ifft(Transfer(10,f,ds,dr,d_x+element_position(i),H,cw,cb,zw,cb),...
                'symmetric');
    ht = ht + ht_temp;
end
plot(conv(ht,CW));
legend("ÇÐÏò","·¨Ïò","¾µÏñ");