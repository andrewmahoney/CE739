%11/17/21
%Lax-Wendroff scheme

clear all
close all
clc

u = 0.08; %set x velocity
dx = 0.008; %x spacing
dt = 0.1; %time increment
y = u*dt/dx %calculate and display gamma
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

%% Lax-Wendroff scheme
for i = 1:length(x)
    if x(i) <= 0.2
        T_lw(1,i) = 1 - (10*x(i) - 1)^2;
    else
        T_lw(1,i) = 0;
    end
end
T_lw(1:length(t),1) = 0; %left end B.C.
T_lw(1:length(t),length(x)) = 0; %right end B.C.

%calculate values, with each row as a time step, each column as a value of x
for i = 2:length(t)
    for j = 2:length(x) - 1
        T_lw(i,j) = T_lw(i-1,j) - y/2*(T_lw(i-1,j+1) - T_lw(i-1,j-1)) + y^2/2*(T_lw(i-1,j+1) - 2*T_lw(i-1,j) + T_lw(i-1,j-1));
    end
end

figure(1)
plot(x_exact,T_exact(find(t_exact==0),:), 'r'), xlabel('x'), ylabel('T'), grid
hold on
plot(x_exact,T_exact(find(t_exact==4),:), 'b')
plot(x_exact,T_exact(find(t_exact==8),:), 'g')

plot(x,T_lw(find(t==0),:), 'r o')
plot(x,T_lw(find(t==4),:), 'b +')
plot(x,T_lw(find(t==8),:), 'g x'), legend('exact: t = 0', 'exact: t = 4', 'exact: t = 8', 'LW: t = 0','LW: t = 4', 'LW: t = 8')
        
        
        
        
        