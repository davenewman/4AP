clear;
clc;
close all;

x = linspace(0,1,500);
U0 = 1;
UL = 0;
A = 1;
k = 10;
L = 1;
v = 1;

%Assignment's soln
u_1_prof = (1 - ( sinh(k*(L-x)) + sinh(k*x))/sinh(k*L))*A/(k^2) + U0*sinh(k*(L-x))/sinh(k*L);

%Assignment's soln
u_2_prof = (1 - (cosh(k*x)/cosh(k*L)))*(A/k^2) - (v/k)*sinh(k*(L-x))/sinh(k*L);

%Analytical Solutions for problem 1
u1_10 = (exp(10*x)*(exp(10) - 101))/(100*(exp(20) - 1)) - (exp(-10*x)*(exp(10) - 101*exp(20)))/(100*(exp(20) - 1)) - 1/100;
u1_100 = (exp(100*x)*(exp(100) - 10001))/(10000*(exp(200) - 1)) - (exp(-100*x)*(exp(100) - 10001*exp(200)))/(10000*(exp(200) - 1)) - 1/10000;

%Analytical Solutions for problem 2
u2_10 = (exp(-10*x)*(exp(10) - 10*exp(20)))/(100*(exp(20) + 1)) + (exp(10*x)*(exp(10) + 10))/(100*(exp(20) + 1)) - 1/100;
u2_100 = (exp(-100*x)*(exp(100) - 100*exp(200)))/(10000*(exp(200) + 1)) + (exp(100*x)*(exp(100) + 100))/(10000*(exp(200) + 1)) - 1/10000;

%comparison of Assignment soln and MATLAB 
plot(x,u_1_prof,x,u1_10);
grid on;
legend('Given in the Assignment','Solution using MATLAB Solver','Location','East');
title('Analytical Solution for k = 10, u(0) = 1, u(1) = 0');

figure;
plot(x,u_2_prof,x,u2_10);
legend('Given in the Assignment','Solution using MATLAB Solver','Location','East');
title("Analytical Solution for k = 10, u'(0) = 1, u(1) = 0");
grid on;


%Problem 1, k = 10
n_10_1 = dlmread('t_0_n_10_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_20_1 = dlmread('t_0_n_20_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_40_1 = dlmread('t_0_n_40_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_80_1 = dlmread('t_0_n_80_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_160_1 = dlmread('t_0_n_160_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_320_1 = dlmread('t_0_n_320_k2_100.0_A_1.0_L_1.0_R_0.0.txt');

%Problem 1, k = 100
n_10_2 = dlmread('t_0_n_10_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_20_2 = dlmread('t_0_n_20_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_40_2 = dlmread('t_0_n_40_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_80_2 = dlmread('t_0_n_80_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_160_2 = dlmread('t_0_n_160_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_320_2 = dlmread('t_0_n_320_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');

%Problem 2, k = 10
n_10_3 = dlmread('t_1_n_10_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_20_3 = dlmread('t_1_n_20_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_40_3 = dlmread('t_1_n_40_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_80_3 = dlmread('t_1_n_80_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_160_3 = dlmread('t_1_n_160_k2_100.0_A_1.0_L_1.0_R_0.0.txt');
n_320_3 = dlmread('t_1_n_320_k2_100.0_A_1.0_L_1.0_R_0.0.txt');

%Problem 2, k = 100
n_10_4 = dlmread('t_1_n_10_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_20_4 = dlmread('t_1_n_20_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_40_4 = dlmread('t_1_n_40_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_80_4 = dlmread('t_1_n_80_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_160_4 = dlmread('t_1_n_160_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');
n_320_4 = dlmread('t_1_n_320_k2_10000.0_A_1.0_L_1.0_R_0.0.txt');


%Figure for problem 1, k = 10
figure;
grid on;
hold on;
plot(x,u1_10,'k',n_10_1(:,1),n_10_1(:,2),'r*');
plot(n_20_1(:,1),n_20_1(:,2),'b*');
plot(n_40_1(:,1),n_40_1(:,2),'m*');
plot(n_80_1(:,1),n_80_1(:,2),'g*');
legend('Analytical Solution','n = 10','n = 20','n = 40','n = 80');
title('Convergence of Approximation for k = 10, u(0) = 1, u(1) = 0');

%Figure for problem 1, k = 100
figure;
grid on;
hold on;
plot(x,u1_100,'k',n_10_2(:,1),n_10_2(:,2),'r*');
plot(n_20_2(:,1),n_20_2(:,2),'b*');
plot(n_40_2(:,1),n_40_2(:,2),'m*');
plot(n_80_2(:,1),n_80_2(:,2),'g*');
legend('Analytical Solution','n = 10','n = 20','n = 40','n = 80');
title('Convergence of Approximation for k = 100, u(0) = 1, u(1) = 0');

%Figure for problem 2, k = 10
figure;
grid on;
hold on;
plot(x,u2_10,'k',n_10_3(:,1),n_10_3(:,2),'r*');
plot(n_20_3(:,1),n_20_3(:,2),'b*');
plot(n_40_3(:,1),n_40_3(:,2),'m*');
plot(n_80_3(:,1),n_80_3(:,2),'g*');
legend('Analytical Solution','n = 10','n = 20','n = 40','n = 80','Location','East');
title("Convergence of Approximation for k = 10, u'(0) = 1, u(1) = 0");

%Figure for problem 2, k = 100
figure;
grid on;
hold on;
plot(x,u2_100,'k',n_10_4(:,1),n_10_4(:,2),'r*');
plot(n_20_4(:,1),n_20_4(:,2),'b*');
plot(n_40_4(:,1),n_40_4(:,2),'m*');
plot(n_80_4(:,1),n_80_4(:,2),'g*');
legend('Analytical Solution','n = 10','n = 20','n = 40','n = 80','Location','East');
title("Convergence of Approximation for k = 100, u'(0) = 1, u(1) = 0");

%Determination of order of accuracy
[order_1_10, differences_1_10] = OrderofAccuracy(n_40_1(:,2),n_80_1(:,2),1);
[order_1_100, differences_1_100] = OrderofAccuracy(n_40_2(:,2),n_80_2(:,2),2);
[order_2_10, differences_2_10] = OrderofAccuracy(n_40_3(:,2),n_80_3(:,2),3);
[order_2_100, differences_2_100] = OrderofAccuracy(n_40_4(:,2),n_80_4(:,2),4);

fprintf('Order of Accuracy:\n');
fprintf('Problem 1: n = 10, Order = %f\n',order_1_10);
fprintf('Problem 1: n = 100, Order = %f\n',order_1_100);
fprintf('Problem 2: n = 10, Order = %f\n',order_2_10);
fprintf('Problem 2: n = 100, Order = %f\n',order_2_100);

figure;
plot(differences_1_10(:,1),differences_1_10(:,2));
title("Error as Function of x for k = 10, u(0) = 1, u(1) = 0");

figure;
plot(differences_1_100(:,1),differences_1_100(:,2));
title("Error as Function of x for k = 100, u(0) = 1, u(1) = 0");

figure;
plot(differences_2_10(:,1),differences_2_10(:,2));
title("Error as Function of x for k = 10, u'(0) = 1, u(1) = 0");

figure;
plot(differences_2_100(:,1),differences_2_100(:,2));
title("Error as Function of x for k = 100, u'(0) = 1, u(1) = 0");

