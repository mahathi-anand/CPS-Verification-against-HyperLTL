%Verifying initial state opacity for vehicle dynamics 
%1st hypersafety element: (a_0,~a_0) to (o,~o) \cup (~o,o)
clc;
clear;
close all;

%============== Parameters for SOS program ===============
%x1- Absolute position of original system 
%v1- Absolute velocity of original system
%x2- Absolute position of virtual system
%v2- Absolute velocity of virtual system 
%u1- control input of original system (assosiciated with universal
%quantifier)
%u2- control input of virtual system (existential quantifier)
tic

epsilon=0;
delta=0.15;

syms x1 x2 v1 v2 u1
vars=[x1;x2;v1;v2;u1];
varx=[x1;x2;v1;v2];

%============= Dynamics of the original system ===========
f_x1= x1+v1+0.5*u1;
f_v1=v1+u1;

%============ Dynamics of the additional system ==========
u2=0.983*v1-v2+u1; %controller

f_x2=x2+v2+0.5*u2;
f_v2=v2+u2;

%=========== Regions of interest for augmented system ==============

%state space

x1_min=0;
x1_max=8;

v1_min=0;
v1_max=0.6; 

x2_min=0;
x2_max=8;

v2_min=0;
v2_max=0.6;

%secret states
%Only for absolute position of the original system.
x10_min=0;
x10_max=1;

%non secret states
%For absolute position of virtual system
x20_min=1;
x20_max=8;

%input constraint 
u_min=-0.04;
u_max=0.04;

% Semialgebraic sets for state space
gx1=(x1-x1_min)*(x1_max-x1);
gx2=(x2-x2_min)*(x2_max-x2);
gv1=(v1-v1_min)*(v1_max-v1);
gv2=(v2-v2_min)*(v2_max-v2);

% Semialgebraic sets for initial condition

g01=(x1-x10_min)*(x10_max-x1); %must belong in secret set
g02=(x2-x20_min)*(x20_max-x2); %must belong in non-secret set
g03=-(v1-v2)^2+delta^2; %outputs must be delta close

g0=[g01; g02; gv1; gv2; g03];

% Semialgebraic sets for unsafe set (the region where outputs don't match)

gun=(v1-v2)^2-delta^2-0.01; 

gu=[gx1; gx2; gv1; gv2; gun];

% Semialgebraic set for input

gi=(u1-u_min)*(u_max-u1);

%the whole augmented state-input space

g=[gx1; gx2; gv1; gv2; gi];

% ========================= Initialization the sum of squares program =========================
prog = sosprogram(vars);

%======== Barrier Certificate and Lagrangian Multipliers ======

%defining monomials 

deg=2;
u_deg=2;
Mon_Barrier= monomials(varx,[0:deg]);

Mon_x1=monomials(x1,[0:deg]);
Mon_x2=monomials(x2,[0:deg]);
Mon_v1=monomials(v1,[0:deg]);
Mon_v2=monomials(v2,[0:deg]);
Mon_vv=monomials([v1;v2],[0:deg]);

Mon_u1=monomials(u1,[0:u_deg]);

%Langrangian for 1st condition (initial set)

[prog,L01]=sospolyvar(prog,Mon_x1,'wscoeff');         
[prog,L02]=sospolyvar(prog,Mon_x2,'wscoeff');    
[prog,L03]=sospolyvar(prog,Mon_v1,'wscoeff');    
[prog,L04]=sospolyvar(prog,Mon_v2,'wscoeff'); 
[prog,L05]=sospolyvar(prog,Mon_vv,'wscoeff');


L0=[L01 L02 L03 L04 L05];

%Langrangian for 2nd condition (unsafe/observable set)

[prog,Lu1]=sospolyvar(prog,Mon_x1,'wscoeff');     
[prog,Lu2]=sospolyvar(prog,Mon_x2,'wscoeff');    
[prog,Lu3]=sospolyvar(prog,Mon_v1,'wscoeff');     
[prog,Lu4]=sospolyvar(prog,Mon_v2,'wscoeff');  
[prog,Lu5]=sospolyvar(prog,Mon_vv,'wscoeff');

Lu=[Lu1 Lu2 Lu3 Lu4 Lu5];

%Lagrangian for last condition (state space and input space)

[prog,L1]=sospolyvar(prog,Mon_x1,'wscoeff'); 
[prog,L2]=sospolyvar(prog,Mon_x2,'wscoeff'); 
[prog,L3]=sospolyvar(prog,Mon_v1,'wscoeff'); 
[prog,L4]=sospolyvar(prog,Mon_v2,'wscoeff'); 
[prog,Li]=sospolyvar(prog,Mon_u1,'wscoeff'); 

L=[L1 L2 L3 L4 Li];

%defining barrier

[prog,Barrier] = sospolyvar(prog,Mon_Barrier,'wscoeff');

B_f = subs(Barrier,{x1,x2,v1,v2},{f_x1,f_x2,f_v1,f_v2});



%========== SOS Constraints ================

%  B<=0 (initial)
prog = sosineq(prog,-Barrier+epsilon-L0*g0); 

%   B>0 (unsafe)
prog = sosineq(prog,Barrier-Lu*gu-0.015-epsilon);

%Bf-B <=0 (barrier decreasing over state space)
prog = sosineq(prog, -B_f+Barrier-L*g);  


%Positiveness of Lagrangian coefficients 

prog = sosineq(prog, L01);  
prog = sosineq(prog, L02); 
prog = sosineq(prog, L03);  
prog = sosineq(prog, L04); 
prog = sosineq(prog, L05);

prog = sosineq(prog, Lu1); 
prog = sosineq(prog, Lu2); 
prog = sosineq(prog, Lu3); 
prog = sosineq(prog, Lu4); 
prog = sosineq(prog, Lu5); 


prog = sosineq(prog, L1); 
prog = sosineq(prog, L2); 
prog = sosineq(prog, L3); 
prog = sosineq(prog, L4); 
prog = sosineq(prog, Li); 

% ========================= Call sos solver (requires SeDuMi installed) =========================
prog = sossolve(prog); 

% ========================= Cheking the results =========================
barrier  = sosgetsol(prog,Barrier)

SOLV2 = sosgetsol(prog,-Barrier-L0*g0+epsilon);      
SOLV3 = sosgetsol(prog, Barrier-Lu*gu-epsilon-0.01);     
SOLV4 = sosgetsol(prog, -B_f+Barrier-L*g);

%========  checking SOS ======== 
[P_b2,z_b2]= findsos(SOLV2);
[P_b3,z_b3]= findsos(SOLV3);
[P_b4,z_b4]= findsos(SOLV4);

aaa = [length(P_b2),length(P_b3),length(P_b4)];
mmm = [length(z_b2),length(z_b3),length(z_b4)]

toc