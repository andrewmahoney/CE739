%% 10/20/21
%problem is y' = e^(y-t), y(0) = y0
%exact solution is y(t) = -ln(e^(-y0) + e^(-t) - 1)

clear all
clc
close all
y0_1 = -1e-5;
y0_2 = -1;
h = 0.2;
t_exact = 0:h:1;
y_exact_1 = -log(exp(-y0_1) + exp(-t_exact) - 1);
y_exact_2 = -log(exp(-y0_2) + exp(-t_exact) - 1);



%% Linearized implicit Euler
%scheme is y(n+1) = (y(n) + (1 + y(n))*h*e^(y(n) - t(n+1)))/(1 - h*e^(y(n) - t(n+1)))

y_lie_1(1) = y0_1; %initial condition of y0 = -1e-5

for j = 2:length(t_exact);
    y_lie_1(j) = (y_lie_1(j-1) + (1 - y_lie_1(j-1))*h*exp(y_lie_1(j-1) - t_exact(j)))/(1 - h*exp(y_lie_1(j-1) - t_exact(j)));
end

error_lie_1 = abs(y_exact_1 - y_lie_1);

figure(1)
plot(t_exact, error_lie_1, 'r +')
hold on
%-------------------------------------------------
y_lie_2(1) = y0_2; %initial condition of y0 = -1

for j = 2:length(t_exact);
    y_lie_2(j) = (y_lie_2(j-1) + (1 - y_lie_2(j-1))*h*exp(y_lie_2(j-1) - t_exact(j)))/(1 - h*exp(y_lie_2(j-1) - t_exact(j)));
end

error_lie_2 = abs(y_exact_2 - y_lie_2);

figure(1)
plot(t_exact, error_lie_2, 'b x'), legend('y0 = -1e-5', 'y0 = -1')

%% Fully implicit Euler
%need to find roots of (y(n+1) - y(n))/h - e^(y(n+1) - t(n+1)) = 0

y_ie_1(1) = y0_1; %initial condition of y0 = -1e-5

%from Excel, we have

y_ie_1 = [y0_1, 0.19999, 0.39999, 0.59999, 0.79999, 0.99999];

%myfun = @(yn1, yn, tn1, t) (yn1 - yn)/t - exp(yn1 - tn1);
%t = 0.2;
%yn = y0;
%tn1 = t_exact(2);
%g = @(yn1) myfun(yn1, yn, tn1, t);
%yn1 = fzero(g,1)

%end


%for j = 2:length(t_exact);
    %y_ie(j) = fzero(@(y_ie1) impe(y_ie1, y_ie(j-1), t_exact(j), 0.2),1);
%end


figure(2)

plot(t_exact,y_exact_1)
hold on
plot(t_exact, y_lie_1, 'r +')
plot(t_exact, y_ie_1, 'b *'), legend('exact solution', 'linearized implicit Euler','fully implicit Euler')

%------------------------------
y_ie_2(1) = y0_2; %initial condition of y0 = -1

%from Excel, we have

y_ie_2 = [y0_2, -0.9356, -0.87999, -0.8322, -0.7915, -0.757];

figure(3)
 
plot(t_exact,y_exact_2)
hold on
plot(t_exact, y_lie_2, 'r +')
plot(t_exact, y_ie_2, 'b *'), legend('exact solution','linearized implicit Euler','fully implicit Euler')