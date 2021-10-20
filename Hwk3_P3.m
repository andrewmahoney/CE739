%% 10/20/21
%problem is y' = iy, y(0) = 1
%f = iy where i = sqrt(-1)
%Exact solution is y = e^(i*t)

clear all
close all

%% Explicit Euler
%scheme is y(n+1) = y(n) + h*i*y(n)

h = 0.1;
t = 0:h:20;
y_exact = exp(i*t); %exact solution
plot(t,y_exact)
hold on

y_ee(1) = 1;

for j = 2:length(t);
    y_ee(j) = y_ee(j-1) + h*i*y_ee(j-1);
end

plot(t,y_ee, 'r +')


%% Runge-Kutta second order
%scheme is ytilde = y(n) + h*i*y(n), y(n+1) = y(n) + h/2*i*y(n) + h/2*i*ytilde

y_rk2(1) = 1;
ytilde(1) = y_rk2(1) + h*i*y_rk2(1);

for j = 2:length(t);
    ytilde(j) = y_rk2(j-1) + h*i*y_rk2(j-1);
    y_rk2(j) = y_rk2(j-1) + h/2*i*y_rk2(j-1) + h/2*i*ytilde(j);
end

plot(t,y_rk2, 'b x')


%% Runge-Kutta fourth order
%see hwk for scheme, cumbersome to type here

y_rk4(1) = 1;

for j = 2:length(t);
    K1(j-1) = h*i*y_rk4(j-1);
    K2(j-1) = h*i*(y_rk4(j-1) + 1/2*K1(j-1));
    K3(j-1) = h*i*(y_rk4(j-1) + 1/2*K2(j-1));
    K4(j-1) = h*i*(y_rk4(j-1) + K3(j-1));
    
    y_rk4(j) = y_rk4(j-1) + 1/6*K1(j-1) + 1/3*(K2(j-1) + K3(j-1)) + 1/6*K4(j-1);
end

plot(t,y_rk4, 'g .')
    

%% Leapfrog Method
%scheme is y(n+1) = y(n-1) + 2*h*i*y(n)

y_lf(1) = 1;
y_lf_tilde(2) = y_lf(1) + h*i*y_lf(1); %use RK2 for 2nd point, since leapfrog method requires 3
y_lf(2) = y_lf(1) + h/2*i*y_lf(1) + h/2*i*y_lf_tilde(2);


for j = 3:length(t);
    y_lf(j) = y_lf(j-2) + 2*h*i*y_lf(j-1);
end



plot(t, y_lf, 'm *')



%% Adams Bashforth (2nd order)
%scheme is y(n+1) = y(n) + 3*h/2*i*y(n) - h/2*i*y(n-1)

y_ab(1) = 1;
y_ab_tilde(2) = y_ab(1) + h*i*y_ab(1); %use RK2 for 2nd point, since Adams Bashforth requires 3
y_ab(2) = y_ab(1) + h/2*i*y_ab(1) + h/2*i*y_ab_tilde(2);

for j = 3:length(t);
    y_ab(j) = y_ab(j-1) + 3*h/2*i*y_ab(j-1) - h/2*i*y_ab(j-2);
end

plot(t,y_ab, 'k s'), legend('exact solution','explicit Euler',...
    'Runge-Kutta 2nd order','Runge-Kutta 4th order','leapfrog','Adams Bashforth')
    
    
    
    
    
    
    