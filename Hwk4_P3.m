%11/16/21

clear all
close all
clc

u = 0.08; %set x velocity
dx = 0.005; %x spacing
dt = 0.1; %time increment
CFL = u*dt/dx %calculate and display CFL number
x = 0:dx:1; %set x range
t = 0:dt:8; %set time range



%% Exact Solution, constant u

dx_exact = 0.01; %x spacing
dt_exact = 0.01; %time increment
x_exact = 0:dx_exact:1; %set x range
t_exact = 0:dt_exact:8; %set time range

%setting initial and boundary conditions for T
for i = 1:length(t_exact)
    for j = 1:length(x_exact)
        if x_exact(j) - u*t_exact(i) >= 0 & x_exact(j) - u*t_exact(i) <= 0.2
            T_exact(i,j) = 1 - (10*(x_exact(j) - u*t_exact(i)) - 1)^2;
        else
            T_exact(i,j) = 0;
        end
    end
end

%% Explicit Euler, constant u

for i = 1:length(x)
    if x(i) <= 0.2
        T_EE(1,i) = 1 - (10*x(i) - 1)^2;
    else
        T_EE(1,i) = 0;
    end
end
T_EE(1:length(t),1) = 0; %left end B.C.
T_EE(1:length(t),length(x)) = 0; %right end B.C.

%calculate values, with each row as a time step, each column as a value of x
for i = 2:length(t)
    for j = 2:length(x) - 1
        T_EE(i,j) = T_EE(i-1,j) - dt*u*(T_EE(i-1,j+1) - T_EE(i-1,j-1))/(2*dx);
    end
end

figure(1)
plot(x_exact,T_exact(find(t_exact==0),:), 'r'), xlabel('x'), ylabel('T'), grid
hold on
plot(x_exact,T_exact(find(t_exact==4),:), 'b')
plot(x_exact,T_exact(find(t_exact==8),:), 'g')

plot(x,T_EE(find(t==0),:), 'r o')
plot(x,T_EE(find(t==4),:), 'b +')
plot(x,T_EE(find(t==8),:), 'g x'), legend('exact: t = 0', 'exact: t = 4', 'exact: t = 8', 'EE: t = 0','EE: t = 4', 'EE: t = 8')


%% Leapfrog, constant u
%setting initial and boundary conditions for T
for i = 1:length(x)
    if x(i) <= 0.2
        T_lf(1,i) = 1 - (10*x(i) - 1)^2;
    else
        T_lf(1,i) = 0;
    end
end
T_lf(1:length(t),1) = 0; %left end B.C.
T_lf(1:length(t),length(x)) = 0; %right end B.C.

%calculate values, with each row as a time step, each column as a value of x
for i = 3:length(t)
    for j = 2:length(x) - 1
        T_lf(i,j) = T_lf(i-2,j) - 2*dt*u*(T_lf(i-1,j+1) - T_lf(i-1,j-1))/(2*dx);
    end
end

figure(2)
plot(x_exact,T_exact(find(t_exact==0),:), 'r'), xlabel('x'), ylabel('T'), grid
hold on
plot(x_exact,T_exact(find(t_exact==4),:), 'b')
plot(x_exact,T_exact(find(t_exact==8),:), 'g')

plot(x,T_lf(find(t==0),:), 'r o')
plot(x,T_lf(find(t==4),:), 'b +')
plot(x,T_lf(find(t==8),:), 'g x'), legend('exact: t = 0', 'exact: t = 4', 'exact: t = 8', 'LF: t = 0','LF: t = 4', 'LF: t = 8')

%% Exact Solution, variable u

dx_exact2 = 0.001; %x spacing
dt_exact2 = 0.001; %time increment
x_exact2 = 0:dx_exact2:1; %set x range
t_exact2 = 0:dt_exact2:8; %set time range

%setting initial and boundary conditions for T
for i = 1:length(t_exact2)
    for j = 1:length(x_exact2)
        u2(j) = 0.2*sin(pi*x_exact2(j));
        if x_exact2(j) - u2(j)*t_exact2(i) >= 0 & x_exact2(j) - u2(j)*t_exact2(i) <= 0.2
            T_exact2(i,j) = 1 - (10*(x_exact2(j) - u2(j)*t_exact2(i)) - 1)^2;
        else
            T_exact2(i,j) = 0;
        end
    end
end

%% Explicit Euler, variable u

dx2 = 0.01; %x spacing
dt2 = 0.01; %time increment
x2 = 0:dx2:1; %set x range
t2 = 0:dt2:8; %set time range

%setting initial and boundary conditions for T
for i = 1:length(x2)
    if x2(i) <= 0.2
        T_EE2(1,i) = 1 - (10*x2(i) - 1)^2;
    else
        T_EE2(1,i) = 0;
    end
end
T_EE2(1:length(t2),1) = 0; %left end B.C.
T_EE2(1:length(t2),length(x2)) = 0; %right end B.C.

%calculate values, with each row as a time step, each column as a value of x
for i = 2:length(t2)
    for j = 2:length(x2) - 1
        u2(j) = 0.2*sin(pi*x2(j));
        T_EE2(i,j) = T_EE2(i-1,j) - dt2*u2(j)*(T_EE2(i-1,j+1) - T_EE2(i-1,j-1))/(2*dx2);
    end
end

figure(3)
plot(x_exact2,T_exact2(find(t_exact2==0),:), 'r'), xlabel('x'), ylabel('T'), grid
hold on
plot(x_exact2,T_exact2(find(t_exact2==4),:), 'b')
plot(x_exact2,T_exact2(find(t_exact2==8),:), 'g')

plot(x2,T_EE2(find(t2==0),:), 'r o')
plot(x2,T_EE2(find(t2==4),:), 'b +')
plot(x2,T_EE2(find(t2==8),:), 'g x'), legend('exact: t = 0', 'exact: t = 4', 'exact: t = 8', 'EE: t = 0','EE: t = 4', 'EE: t = 8')
%setting initial and boundary conditions for T

%% Leapfrog, variable u
%setting initial and boundary conditions for T
for i = 1:length(x2)
    if x2(i) <= 0.2
        T_lf2(1,i) = 1 - (10*x2(i) - 1)^2;
    else
        T_lf2(1,i) = 0;
    end
end
T_lf2(1:length(t2),1) = 0; %left end B.C.
T_lf2(1:length(t2),length(x2)) = 0; %right end B.C.

%calculate values, with each row as a time step, each column as a value of x
for i = 3:length(t2)
    for j = 2:length(x2) - 1
        u2(j) = 0.2*sin(pi*x2(j));
        T_lf2(i,j) = T_lf2(i-2,j) - 2*dt2*u2(j)*(T_lf2(i-1,j+1) - T_lf2(i-1,j-1))/(2*dx2);
    end
end

figure(4)
plot(x_exact2,T_exact2(find(t_exact2==0),:), 'r'), xlabel('x'), ylabel('T'), grid
hold on
plot(x_exact2,T_exact2(find(t_exact2==4),:), 'b')
plot(x_exact2,T_exact2(find(t_exact2==8),:), 'g')

plot(x2,T_lf2(find(t2==0),:), 'r o')
plot(x2,T_lf2(find(t2==4),:), 'b +')
plot(x2,T_lf2(find(t2==8),:), 'g x'), legend('exact: t = 0', 'exact: t = 4', 'exact: t = 8', 'LF: t = 0','LF: t = 4', 'LF: t = 8')

