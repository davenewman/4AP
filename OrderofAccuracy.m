function [order, differences] = OrderofAccuracy(n_approximation,n2_approximation, flag)

x = linspace(0,1,length(n_approximation))';
x2 = linspace(0,1,length(n2_approximation))';
differences = zeros(length(x2),2);

if flag == 1
    u1 = (exp(10*x)*(exp(10) - 101))/(100*(exp(20) - 1)) - (exp(-10*x)*(exp(10) - 101*exp(20)))/(100*(exp(20) - 1)) - 1/100;
    u2 = (exp(10*x2)*(exp(10) - 101))/(100*(exp(20) - 1)) - (exp(-10*x2)*(exp(10) - 101*exp(20)))/(100*(exp(20) - 1)) - 1/100;
elseif flag == 2
    u1 = (exp(100*x)*(exp(100) - 10001))/(10000*(exp(200) - 1)) - (exp(-100*x)*(exp(100) - 10001*exp(200)))/(10000*(exp(200) - 1)) - 1/10000;
    u2 = (exp(100*x2)*(exp(100) - 10001))/(10000*(exp(200) - 1)) - (exp(-100*x2)*(exp(100) - 10001*exp(200)))/(10000*(exp(200) - 1)) - 1/10000;
elseif flag == 3
    u1 = (exp(-10*x)*(exp(10) - 10*exp(20)))/(100*(exp(20) + 1)) + (exp(10*x)*(exp(10) + 10))/(100*(exp(20) + 1)) - 1/100;
    u2 = (exp(-10*x2)*(exp(10) - 10*exp(20)))/(100*(exp(20) + 1)) + (exp(10*x2)*(exp(10) + 10))/(100*(exp(20) + 1)) - 1/100;
else
    u1 = (exp(-100*x)*(exp(100) - 100*exp(200)))/(10000*(exp(200) + 1)) + (exp(100*x)*(exp(100) + 100))/(10000*(exp(200) + 1)) - 1/10000;
    u2 = (exp(-100*x2)*(exp(100) - 100*exp(200)))/(10000*(exp(200) + 1)) + (exp(100*x2)*(exp(100) + 100))/(10000*(exp(200) + 1)) - 1/10000;
end



order = log( max( abs(n_approximation - u1) )/ max( abs(n2_approximation - u2) ) )/log(2);

%only looking at differences in the solution with 80 points
differences(:,1) = x2';
differences(:,2) = abs(n2_approximation - u2);