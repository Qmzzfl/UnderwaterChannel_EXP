clc;close all;clear;
t_end = 10;
dt = 0.001;
t = 0:dt:t_end;
f = 200;
me = 1;

measure_number = 500;

measure_signal = zeros(measure_number,length(t));
df_range = 10;
ddf = df_range*2/measure_number;
i = 1;
for df = -df_range:ddf:df_range-ddf
    measure_signal(i,:) = sin((f+me*t+df).*t);
    i = i + 1;
end

t_r_end = 20;
t_r = 0:dt:t_r_end;
df_r = 10;
received_signal = sin((f+me*t_r+df_r).*t_r);
measure_result = zeros(measure_number,length(t)+length(t_r)-1);
for j = 1:measure_number
    measure_result(j,:) = xcorr2(received_signal,measure_signal(j,:));
end

surf(measure_result);
shading interp