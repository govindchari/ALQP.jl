clear;clc;close all
n_list = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
osqp = [3.93125e-6, 4.850285714285715e-6, 7.317e-6, 1.342e-5, 3.385e-5, 0.00010952, 0.000394529, 0.001499322, 0.005815006, 0.030798397];
alqp = [1.6763866666666667e-7, 3.712621359223301e-7, 1.1405e-6, 8.00175e-6, 4.6639e-5, 0.000221191, 0.00096504, 0.004001479, 0.017896095, 0.073445111];
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