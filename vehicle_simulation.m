% Manual Computation of Controller

%state space
  clear all
  clc

x1_min=0;
x1_max=8;

x2_min=0;
x2_max=8;

v1_min=0;
v1_max=0.6;

v2_min=0;
v2_max=0.6;


%initial (secret) states
x10_min=0;
x10_max=1;

v10_min=0;
v10_max=0.6;

x20_min=1;
x20_max=8;

v20_min=0;
v20_max=0.6;


%input 

u_min=-0.04;
u_max=0.04;

delta=0.15;

c = get (0, 'DefaultAxesColorOrder' );
for p=1:50

%Generate random initial state 

x1_0 = (x10_max-x10_min).*rand(1)+x10_min;
x2_0 = (x20_max-x20_min).*rand(1)+x20_min;

v1_0 = (v10_max-v10_min).*rand(1)+v10_min;
v2_0 = (v10_max-v10_min).*rand(1)+v10_min;

while (v1_0-v2_0)^2>delta^2
v1_0 = (v10_max-v10_min).*rand(1)+v10_min;
v2_0 = (v10_max-v10_min).*rand(1)+v10_min;
end

L = 50;

x1=zeros(1,L);
x2=zeros(1,L);
v1=zeros(1,L);
v2=zeros(1,L);

x1(1)=x1_0;
x2(1)=x2_0;
v1(1)=v1_0;
v2(1)=v2_0;
u1=zeros(1,L-1);
%define dynamics


for i=1:L-1
    u1(i)=u_min+rand(1)*(u_max-u_min);
    u2(i)=0.983*v1(i)-v2(i)+u1(i);

    x1(i+1)= x1(i)+v1(i)+0.5*u1(i);
    v1(i+1)=v1(i)+u1(i);
    x2(i+1)=x2(i)+v2(i)+0.5*u2(i);
    v2(i+1)=v2(i)+u2(i);
    
    if x1(i+1) >= x1_max
    x1(i+1) = x1_max;
    end
    if x1(i+1) <= x1_min
        x1(i+1) = x1_min;
    end
    if x2(i+1) >= x2_max
    x2(i+1) =x2_max;
    end
    if x2(i+1) <= x2_min
        x2(i+1) = x2_min;
    end
    if v1(i+1) >= v1_max
    v1(i+1) = v1_max;
    end
    if v1(i+1) <= v1_min
        v1(i+1) = v1_min;
    end
    if v2(i+1) >= v2_max 
    v2(i+1) = v2_max;
    end
    if v2(i+1) <= v2_min
       v2(i+1) = v2_min;
    end
    
end
figure(1);

plot(v1(1:L),v2(1:L),'color',c(mod(p,size(c,1))+1,:));
plot(v1_0, v2_0, '*','color',c(mod(p,size(c,1))+1,:));
xlabel('$v_1$', 'Interpreter', 'latex', 'FontSize',20,'Fontname','Arial');
ylabel('$v_2$', 'Interpreter', 'latex', 'FontSize',20,'Fontname','Arial');
hold on

figure(2);
plot(x1_0, x2_0, '*','color',c(mod(p,size(c,1))+1,:));
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize',20,'Fontname','Arial');
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize',20,'Fontname','Arial');
hold on
end
figure(1);
%plotting unsafe regions
fill( [0.15 0.6 0.6], [0  0.45  0],'c', 'facealpha',0.2,'edgealpha',0);
hold on
fill( [0 0.45 0], [0.15  0.6  0.6],'c', 'facealpha',0.2,'edgealpha',0);
hold off


