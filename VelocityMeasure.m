clc;close all;clear;

% simulation time
t_end = 5;
dt = 0.001;
t = 0:dt:t_end;
conv_t = -t_end:dt:t_end;
fs = 1/dt;

% signal configure
f_start = 200;% static start frequency
f_end = 500;% static  end  frequency
signal = chirp(t,f_start,t(end),f_end);

% Dopler configure
f_shift_lim = 10; % the max of |Dopler frequency shift|
measure_num = 100; % the number of tests
f_shift_axis = linspace(-f_shift_lim,f_shift_lim,measure_num);
f_start_shift = linspace(f_start-f_shift_lim,f_start+f_shift_lim,measure_num);
f_end_shift = linspace(f_end-f_shift_lim,f_end+f_shift_lim,measure_num);

% Dopler signal generate
test_signal = zeros(measure_num,length(t));
for i =1:100
    test_signal(i,:) = chirp(t,f_start_shift(i),t(end),f_end_shift(i));
end

measure_result = zeros(measure_num,2*length(t)-1);
for i = 1:100
   measure_result(i,:) = xcorr(signal,test_signal(i,:));
end

surf(conv_t,f_shift_axis,measure_result);
shading interp