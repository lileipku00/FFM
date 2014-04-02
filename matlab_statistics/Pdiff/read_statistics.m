cd './GSN'
load('t_shift_array_1.mat')
load('t_shift_array_2.mat')
load('t_shift_array_3.mat')
load('t_shift_array_4.mat')
load('t_shift_array_5.mat')
load('t_shift_array_6.mat')
load('t_shift_array_7.mat')
load('t_shift_array_8.mat')

interval = [-5:0.5:5];

t1 = hist(t_shift_array_1, interval);
t2 = hist(t_shift_array_2, interval);
t3 = hist(t_shift_array_3, interval);
t4 = hist(t_shift_array_4, interval);
t5 = hist(t_shift_array_5, interval);
t6 = hist(t_shift_array_6, interval);
t7 = hist(t_shift_array_7, interval);
t8 = hist(t_shift_array_8, interval);
plot(interval, t1)
hold on
plot(interval, t2)
plot(interval, t3)
plot(interval, t4)
plot(interval, t5)
plot(interval, t6)
plot(interval, t7)
plot(interval, t8)

cd '../NO_GSN'

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
plot(interval, t1, 'color', 'red')
plot(interval, t2, 'color', 'red')
plot(interval, t3, 'color', 'red')
plot(interval, t4, 'color', 'red')
plot(interval, t5, 'color', 'red')
plot(interval, t6, 'color', 'red')
plot(interval, t7, 'color', 'red')
plot(interval, t8, 'color', 'red')

