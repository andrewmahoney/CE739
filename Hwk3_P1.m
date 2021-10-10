%% 10/10/21
clear all
%dy/dt = -0.5y
%y(0) = 1 for 0 <= t <=20
%Exact y(t) = exp(-0.5t)

dt = 0.5;
t_exact = 0:dt:20;
Exact_y = exp(-0.5*t_exact);
plot(t_exact,Exact_y);
hold on


%% Explicit Euler
% dt = 1.0
dt = 1;
t = 0:dt:20;
N = length(t);
y_Explicit_1(1) = 1;
%Loop all timesteps
for i = 2:N
    y_Explicit_1(i) = y_Explicit_1(i-1)*(1-0.5*dt);
end
plot(t,y_Explicit_1, 'r o')

% dt  = 4.2
dt = 4.2;
t = 0:dt:20;
N = length(t);
y_Explicit_2(1) = 1;
%Loop all timesteps
for i = 2:N
    y_Explicit_2(i) = y_Explicit_2(i-1)*(1-0.5*dt);
end
plot(t,y_Explicit_2, 'b o')

%% Implicit Euler
dt = 1.0;
t = 0:dt:20;
N = length(t);
y_Implicit_1(1) = 1;
%Loop for all timesteps
for i = 2:N
    y_Implicit_1(i) = y_Implicit_1(i-1)/(1+0.5*dt);
end
plot(t,y_Implicit_1, 'r x');

%dt = 4.2
dt = 4.2;
t = 0:dt:20;
N = length(t);
y_Implicit_2(1) = 1;
%Loop for all timesteps
for i = 2:N
    y_Implicit_2(i) = y_Implicit_2(i-1)/(1+0.5*dt);
end
plot(t,y_Implicit_2, 'b x'), xlim([0 20]), legend("exact", "explicit, h=1.0",...
    "explicit, h=4.2", "implicit, h=1.0", "implicit, h=4.2"), xlabel("t"), ylabel("y(t)")

    
    
    
    
    