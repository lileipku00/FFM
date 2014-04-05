% Simple script to plot the statistical analysis results in MATALB

clear; clc;

% ------------- INPUT --------------------
interval = -6:0.5:6;
% ----------------------------------------

% Following address has the patchline.m which is used here
addpath('~/Work/Scripts/gitHUB/UTILS/OTHERS')

% First plotting GSN ones
cd './GSN'
load('t_shift_array_1.mat')
load('t_shift_array_2.mat')
load('t_shift_array_3.mat')
load('t_shift_array_4.mat')
load('t_shift_array_5.mat')
load('t_shift_array_6.mat')
load('t_shift_array_7.mat')
load('t_shift_array_8.mat')

t1 = hist(t_shift_array_1, interval);
t2 = hist(t_shift_array_2, interval);
t3 = hist(t_shift_array_3, interval);
t4 = hist(t_shift_array_4, interval);
t5 = hist(t_shift_array_5, interval);
t6 = hist(t_shift_array_6, interval);
t7 = hist(t_shift_array_7, interval);
t8 = hist(t_shift_array_8, interval);

patchline(interval, t1, 'edgecolor', [0.478, 0.063, 0.894], ...
    'linewidth', 4, 'edgealpha', 1.0);
hold on
patchline(interval, t2, 'edgecolor', [0., 0., 1.0], ...
    'linewidth', 4, 'edgealpha', 1.0);
patchline(interval, t3, 'edgecolor', [0.043, 0.518, 0.78], ...
    'linewidth', 4, 'edgealpha', 1.0);
patchline(interval, t4, 'edgecolor', [0., 0.498, 0.], ...
    'linewidth', 4, 'edgealpha', 1.0)
patchline(interval, t5, 'edgecolor', [0.749, 0.749, 0.], ...
    'linewidth', 4, 'edgealpha', 1.0)
patchline(interval, t6, 'edgecolor', [0.871, 0.49, 0.894], ...
    'linewidth', 4, 'edgealpha', 1.0)
patchline(interval, t7, 'edgecolor', [0.847, 0.161, 0.], ...
    'linewidth', 4, 'edgealpha', 1.0)
patchline(interval, t8, 'edgecolor', [0.6, 0.2, 0.], ...
    'linewidth', 4, 'edgealpha', 1.0)


p1 = plot(interval, t1, 'color', [0.478, 0.063, 0.894], 'linewidth', 4);
hold on
p2 = plot(interval, t2, 'color', [0., 0., 1.0], 'linewidth', 4);
p3 = plot(interval, t3, 'color', [0.043, 0.518, 0.78], 'linewidth', 4);
p4 = plot(interval, t4, 'color', [0., 0.498, 0.], 'linewidth', 4);
p5 = plot(interval, t5, 'color', [0.749, 0.749, 0.], 'linewidth', 4);
p6 = plot(interval, t6, 'color', [0.871, 0.49, 0.894], 'linewidth', 4);
p7 = plot(interval, t7, 'color', [0.847, 0.161, 0.], 'linewidth', 4);
p8 = plot(interval, t8, 'color', [0.6, 0.2, 0.], 'linewidth', 4);

legend([p1 p2 p3 p4 p5 p6 p7 p8], {'30.0s', '21.2s', '15.0s', '10.6s',...
    '7.5s', '5.3s', '3.7s', '2.7s'})


% Second plotting NO_GSN ones
% cd '../NO_GSN'
% load('t_shift_array_1.mat')
% load('t_shift_array_2.mat')
% load('t_shift_array_3.mat')
% load('t_shift_array_4.mat')
% load('t_shift_array_5.mat')
% load('t_shift_array_6.mat')
% load('t_shift_array_7.mat')
% load('t_shift_array_8.mat')
% 
% t1 = hist(t_shift_array_1, interval);
% t2 = hist(t_shift_array_2, interval);
% t3 = hist(t_shift_array_3, interval);
% t4 = hist(t_shift_array_4, interval);
% t5 = hist(t_shift_array_5, interval);
% t6 = hist(t_shift_array_6, interval);
% t7 = hist(t_shift_array_7, interval);
% t8 = hist(t_shift_array_8, interval);
% patchline(interval, t1, 'linestyle', '--', 'edgecolor', [0.478, 0.063, 0.894], ...
%     'linewidth', 4, 'edgealpha', 0.5)
% hold on
% patchline(interval, t2, 'linestyle', '--', 'edgecolor', [0., 0., 1.0], ...
%     'linewidth', 4, 'edgealpha', 0.5)
% patchline(interval, t3, 'linestyle', '--', 'edgecolor', [0.043, 0.518, 0.78], ...
%     'linewidth', 4, 'edgealpha', 0.5)
% patchline(interval, t4, 'linestyle', '--', 'edgecolor', [0., 0.498, 0.], ...
%     'linewidth', 4, 'linestyle', '--', 'edgealpha', 1.0)
% patchline(interval, t5, 'linestyle', '--', 'edgecolor', [0.749, 0.749, 0.], ...
%     'linewidth', 4, 'edgealpha', 0.5)
% patchline(interval, t6, 'linestyle', '--', 'edgecolor', [0.871, 0.49, 0.894], ...
%     'linewidth', 4, 'edgealpha', 0.5)
% patchline(interval, t7, 'linestyle', '--', 'edgecolor', [0.847, 0.161, 0.], ...
%     'linewidth', 4, 'edgealpha', 0.5)
% patchline(interval, t8, 'linestyle', '--', 'edgecolor', [0.6, 0.2, 0.], ...
%     'linewidth', 4, 'edgealpha', 0.5)

cd '..'