%% 10/10/21
%problem is y'+y(1-y) = 0, y(0) = 1/2
%f = y^2 - y
%Exact solution is y = 1/(1+e^t)

clear all
close all
dt = 0.01;
t_exact = 0:dt:1; %t values for exact solution
Exact_y = 1./(1+exp(t_exact)); %exact solution
figure(1)
plot(t_exact,Exact_y);
hold on

%% Linearized Trapezoidal method ***PLOTTED IN RED***
%Trapezoidal method: y(n+1) = y(n) + h/2*(f(n+1) + f(n))
%Linearization: f(y(n+1),t(n+1) == f(y(n),t(n+1)) + (y(n+1) - y(n))*delf/dely

dt = 1; %t increment
y_approx_a1(1) = 1/2; %initial condition of y(0) = 1/2
t1 = 0:dt:1; %set t values
    for j = 2:length(t1); %loop over values of t
        y_approx_a1(j) = (y_approx_a1(j-1)*(1-dt/2))./(1 - dt*y_approx_a1(j-1) + dt/2);
    end

dt = 0.1;    
y_approx_a2(1) = 1/2;
t2 = 0:dt:1;
    for j = 2:length(t2);
        y_approx_a2(j) = (y_approx_a2(j-1)*(1-dt/2))./(1 - dt*y_approx_a2(j-1) + dt/2);
    end

dt = 0.01;    
y_approx_a3(1) = 1/2;
t3 = 0:dt:1;
    for j = 2:length(t3);
        y_approx_a3(j) = (y_approx_a3(j-1)*(1-dt/2))./(1 - dt*y_approx_a3(j-1) + dt/2);
    end
    
dt = 0.001;    
y_approx_a4(1) = 1/2;
t4 = 0:dt:1;
    for j = 2:length(t4);
        y_approx_a4(j) = (y_approx_a4(j-1)*(1-dt/2))./(1 - dt*y_approx_a4(j-1) + dt/2);
    end
    
    
dt = 0.0001;    
y_approx_a5(1) = 1/2;
t5 = 0:dt:1;
    for j = 2:length(t5);
        y_approx_a5(j) = (y_approx_a5(j-1)*(1-dt/2))./(1 - dt*y_approx_a5(j-1) + dt/2);
    end
    
dt = 0.00001;    
y_approx_a6(1) = 1/2;
t6 = 0:dt:1;
    for j = 2:length(t6);
        y_approx_a6(j) = (y_approx_a6(j-1)*(1-dt/2))./(1 - dt*y_approx_a6(j-1) + dt/2);
    end
    
%plotting
subplot(2,3,1)
plot(t_exact,Exact_y);
hold on
title('h = 1')
plot(t1,y_approx_a1, 'r x')

subplot(2,3,2)
plot(t_exact,Exact_y);
hold on
title('h = 0.1')
plot(t2,y_approx_a2, 'r x')

subplot(2,3,3)
plot(t_exact,Exact_y);
hold on
title('h = 0.01')
plot(t3,y_approx_a3, 'r x')
    
subplot(2,3,4)
plot(t_exact,Exact_y);
hold on
title('h = 0.001')
plot(t4,y_approx_a4, 'r x')

subplot(2,3,5)
plot(t_exact,Exact_y);
hold on
title('h = 0.0001')
plot(t5,y_approx_a5, 'r x')

subplot(2,3,6)
plot(t_exact,Exact_y);
hold on
title('h = 0.00001')
plot(t6,y_approx_a6, 'r x')
    
%% Direct Trapezoidal Method ***PLOTTED IN BLUE***
%Scheme is found d, directly by solving quadratic formula, otherwise same as before
dt = 1;    
y_approx_b1(1) = 1/2;
t1 = 0:dt:1;
    for j = 2:length(t1);
        y_approx_b1(j) = ((1+dt/2) - sqrt((1+dt/2)^2 - 4*(dt/2)*(y_approx_b1(j-1) - dt/2*y_approx_b1(j-1)*(1-y_approx_b1(j-1)))))/dt;
    end
    
dt = 0.1;    
y_approx_b2(1) = 1/2;
t2 = 0:dt:1;
    for j = 2:length(t2);
        y_approx_b2(j) = ((1+dt/2) - sqrt((1+dt/2)^2 - 4*(dt/2)*(y_approx_b2(j-1) - dt/2*y_approx_b2(j-1)*(1-y_approx_b2(j-1)))))/dt;
    end
    
dt = 0.01;    
y_approx_b3(1) = 1/2;
t3 = 0:dt:1;
    for j = 2:length(t3);
        y_approx_b3(j) = ((1+dt/2) - sqrt((1+dt/2)^2 - 4*(dt/2)*(y_approx_b3(j-1) - dt/2*y_approx_b3(j-1)*(1-y_approx_b3(j-1)))))/dt;
    end
    
dt = 0.001;    
y_approx_b4(1) = 1/2;
t4 = 0:dt:1;
    for j = 2:length(t4);
        y_approx_b4(j) = ((1+dt/2) - sqrt((1+dt/2)^2 - 4*(dt/2)*(y_approx_b4(j-1) - dt/2*y_approx_b4(j-1)*(1-y_approx_b4(j-1)))))/dt;
    end
    
dt = 0.0001;    
y_approx_b5(1) = 1/2;
t5 = 0:dt:1;
    for j = 2:length(t5);
        y_approx_b5(j) = ((1+dt/2) - sqrt((1+dt/2)^2 - 4*(dt/2)*(y_approx_b5(j-1) - dt/2*y_approx_b5(j-1)*(1-y_approx_b5(j-1)))))/dt;
    end
    
dt = 0.00001;    
y_approx_b6(1) = 1/2;
t6 = 0:dt:1;
    for j = 2:length(t6);
        y_approx_b6(j) = ((1+dt/2) - sqrt((1+dt/2)^2 - 4*(dt/2)*(y_approx_b6(j-1) - dt/2*y_approx_b6(j-1)*(1-y_approx_b6(j-1)))))/dt;
    end
    
%plotting
figure(2)
subplot(2,3,1)
plot(t_exact,Exact_y);
hold on
title('h = 1')
plot(t1,y_approx_b1, 'b o')

subplot(2,3,2)
plot(t_exact,Exact_y);
hold on
title('h = 0.1')
plot(t2,y_approx_b2, 'b o')

subplot(2,3,3)
plot(t_exact,Exact_y);
hold on
title('h = 0.01')
plot(t3,y_approx_b3, 'b o')

subplot(2,3,4)
plot(t_exact,Exact_y);
hold on
title('h = 0.001')
plot(t4,y_approx_b4, 'b o')

subplot(2,3,5)
plot(t_exact,Exact_y);
hold on
title('h = 0.0001')
plot(t5,y_approx_b5, 'b o')

subplot(2,3,6)
plot(t_exact,Exact_y);
hold on
title('h = 0.00001')
plot(t6,y_approx_b6, 'b o')

figure(3)
plot(t_exact,Exact_y)
hold on
plot(t2, y_approx_a2, 'r x')
plot(t2, y_approx_b2, 'b o')
legend('exact solution','trapezoidal method','direct method')


