clear;clc;close all
n_list = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
osqp = [3.93125e-6, 4.850285714285715e-6, 7.317e-6, 1.342e-5, 3.385e-5, 0.00010952, 0.000394529, 0.001499322, 0.005815006, 0.030798397];
alqp = [1.6426857887874838e-7, 3.711116504854369e-7, 1.1514000000000001e-6, 4.416571428571428e-6, 1.8482e-5, 7.7334e-5, 0.000331087, 0.002462649, 0.019389739, 0.082974987];
naive_alqp = [3.4177e-5, 4.0897e-5, 5.8164e-5, 9.6879e-5, 0.00027465, 0.001395482, 0.007213461, 0.058447229, 0.203651271, 0.75200751];


loglog(n_list, osqp, 'LineWidth',2)
hold on
loglog(n_list, alqp, 'LineWidth',2)
loglog(n_list, naive_alqp, 'LineWidth',2)

xlabel('Number of Assets')
ylabel('Solve Time (sec)')
legend('OSQP', 'ALQP', 'Naive ALQP', 'Location','best')
title('Solve Time vs Number of Assets')
grid on