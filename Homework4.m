clc;clear;close all;

load('E:\ProjectFile\matlab\UnderwaterChannel_EXP\NOF1\mat\NOF1_001.mat')
pcolor(10*log10(abs(h).^2))
shading interp