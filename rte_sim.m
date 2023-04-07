% room temperature regulation of a single room!

%dynamics
%x(k + 1) = x(k) + τs(αe(Te − x(k)) + αH(Th − x(k))u(k))

%parameters
close all
clear all
clc

t_s= 5; %sampling time

T_h= 55; %heater temperature

T_e= 15; %ambient temperature 

a_e= 0.008; %heat exchange coefficients. 

a_h= 0.0036;

l=50;

for p=1:30
x=zeros(1,l);

x(1)= 21+3*rand(1);

u=zeros(1,l);

for i=1:l
    
%u(i)=-1.018e-6*x(i)^4 + 7.563e-5*x(i)^3 - 0.001872*x(i)^2 + 0.02022*x(i) + 0.3944;
u(i)=-0.002398*x(i) + 0.5357;
x(i+1)= x(i)+t_s*a_e*(T_e-x(i)) + a_h*(T_h-x(i))*u(i)*t_s;

end

plot([1:l+1],x);

hold on
end
hold off