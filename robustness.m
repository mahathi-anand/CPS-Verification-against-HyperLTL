%Verifying initial-state robustness 

%quantifiers: \forall \forall 
clc;
clear;
close all;

%============== Parameters for SOS program ===============
%x1- temperature of original system 
%x2-  temperature of virtual system
%u1- control input of original system (assosiciated with universal
%quantifier)
%u2- control input of virtual system (universal quantifier)


%Shifted the dynamics of the system by -20 degrees in order to center the
%system at origin.
tic 

epsilon=0;

syms x1 x2 
vars=[x1;x2];
varx=[x1;x2];

%parameters

t_s= 5;%sampling time

T_h= 55-20; %heater temperature

T_e= 15-20; %ambient temperature 

a_e= 0.008; %heat exchange coefficients. 

a_h= 0.0036;

u1=-0.002398*(x1) + 0.5357;
u2=-0.002398*(x2) + 0.5357;

%============= Dynamics of the original system ===========

f_x1= (x1)+t_s*(a_e*(T_e-(x1))) + a_h*(T_h-(x1))*u1*t_s;

%============ Dynamics of the additional system ==========
f_x2= (x2)+t_s*(a_e*(T_e-(x2))) + a_h*(T_h-(x2))*u2*t_s;

%==== Regions of interest for augmented system ==============

%state space

x1_min=20-20;
x1_max=35-20;

x2_min=20-20;
x2_max=35-20;


%initial states
%Only for temperature of the original system.
x10_min=21-20;
x10_max=21-20;

%initial states
%For temperature of virtual system
x20_min=21-0.5-20;
x20_max=21+0.5-20;

%unsafe states 1

x1u1_min=25-20;
x1u1_max=35-20;

x2u1_min=25-20;
x2u1_max=35-20;

% %unsafe states 2;
% 
 x1u2_min=25-20;
 x1u2_max=35-20;
% 
 x2u2_min=20-20;
 x2u2_max=25-20;
% 
% %unsafe states 3;

x1u3_min=20-20;
x1u3_max=25-20;

x2u3_min=25-20;
x2u3_max=35-20;


% Semialgebraic sets for state space
gx1=(x1-x1_min)*(x1_max-x1);
gx2=(x2-x2_min)*(x2_max-x2);

% Semialgebraic sets for initial condition

g01=(x1-x10_min)*(x10_max-x1); %must belong in initial state
g02=(x2-x20_min)*(x20_max-x2); %must belong in inital set

g0=[g01; g02];

%Semialgebraic sets for unsafe set 1 
gu11=(x1-x1u1_min)*(x1u1_max-x1);
gu12=(x2-x2u1_min)*(x2u1_max-x2);

gu1=[gu11; gu12];

% % Semialgebraic sets for unsafe set 2 
gu21=(x1-x1u2_min)*(x1u2_max-x1);
gu22=(x2-x2u2_min)*(x2u2_max-x2);

gu2=[gu21; gu22];
 
% Semialgebraic sets for unsafe set 3 

gu31=(x1-x1u3_min)*(x1u3_max-x1);
gu32=(x2-x2u3_min)*(x2u3_max-x2);

gu3=[gu31; gu32];


%the whole augmented state-input space

g=[gx1; gx2];

% ========================= Initialization the sum of squares program =========================
prog = sosprogram(vars);

%======== Barrier Certificate and Lagrangian Multipliers ======

%defining monomials 
b_deg=2;
deg=2;

Mon_Barrier= monomials(varx,[0:b_deg]);

Mon_x1=monomials(x1,[0:deg]);
Mon_x2=monomials(x2,[0:deg]);

%Langrangian for 1st condition (initial set)

[prog,L01]=sospolyvar(prog,Mon_x1,'wscoeff');         
[prog,L02]=sospolyvar(prog,Mon_x2,'wscoeff');    

L0=[L01 L02];

%Langrangian for 1st unsafe condition 

 [prog,L03]=sospolyvar(prog,Mon_x1,'wscoeff');    
 [prog,L04]=sospolyvar(prog,Mon_x2,'wscoeff'); 
 
 Lu1=[L03 L04];

% %Langrangian for 2nd unsafe condition 
% 
[prog,L05]=sospolyvar(prog,Mon_x1,'wscoeff');    
[prog,L06]=sospolyvar(prog,Mon_x2,'wscoeff'); 

Lu2=[L05 L06];
% 
% %Langrangian for 3rd unsafe condition 
% 
[prog,L07]=sospolyvar(prog,Mon_x1,'wscoeff');    
[prog,L08]=sospolyvar(prog,Mon_x2,'wscoeff'); 
 
Lu3=[L07 L08];

%Lagrangian for last condition (state space and input space)

[prog,L1]=sospolyvar(prog,Mon_x1,'wscoeff'); 
[prog,L2]=sospolyvar(prog,Mon_x2,'wscoeff'); 


L=[L1 L2];

%defining barrier

[prog,Barrier] = sospolyvar(prog,Mon_Barrier,'wscoeff');

 B_f = subs(Barrier,{x1,x2},{f_x1,f_x2});



%========== SOS Constraints ================
%prog=sosineq(prog,Barrier);

%  B<=0 (initial)
prog = sosineq(prog,-Barrier+epsilon-L0*g0); 

%   B>0 (unsafe)
prog = sosineq(prog,Barrier-Lu1*gu1-0.001-epsilon);
prog = sosineq(prog,Barrier-Lu2*gu2-0.001-epsilon);
prog = sosineq(prog,Barrier-Lu3*gu3-0.001-epsilon);


%Bf-B <=0 (barrier decreasing over state space)
prog = sosineq(prog, -B_f+Barrier-L*g);  


%Positiveness of Lagrangian coefficients 
%initial
prog = sosineq(prog, L01);  
prog = sosineq(prog, L02); 

%unsafe

prog = sosineq(prog, L03);  
prog = sosineq(prog, L04); 
prog = sosineq(prog, L05);
prog = sosineq(prog, L06);
prog = sosineq(prog, L07);
prog = sosineq(prog, L08);

%last condition
prog = sosineq(prog, L1); 
prog = sosineq(prog, L2); 



% ========================= Call sos solver (requires SeDuMi installed) =========================
prog = sossolve(prog); 

% ========================= Cheking the results =========================
barrier  = sosgetsol(prog,Barrier)

SOLV2 = sosgetsol(prog,-Barrier-L0*g0+epsilon);      
SOLV31 = sosgetsol(prog, Barrier-Lu1*gu1-epsilon-0.001);  
SOLV32 = sosgetsol(prog, Barrier-Lu2*gu2-epsilon-0.001);  
SOLV33 = sosgetsol(prog, Barrier-Lu3*gu3-epsilon-0.001);  
% 
SOLV4 = sosgetsol(prog, -B_f+Barrier-L*g);

%========  checking SOS ======== 
[P_b2,z_b2]= findsos(SOLV2);
[P_b31,z_b31]= findsos(SOLV31);
[P_b32,z_b32]= findsos(SOLV32);
[P_b33,z_b33]= findsos(SOLV33);

[P_b4,z_b4]= findsos(SOLV4);

aaa = [length(P_b2),length(P_b31),length(P_b32),length(P_b33),length(P_b4)];
mmm = [length(z_b2),length(z_b31),length(z_b32),length(z_b33),length(z_b4)]

lastcond=sosgetsol(prog,-B_f+Barrier);

%Barrier certificate for unshifted system 

barr=subs(barrier,{x1,x2},{x1-20,x2-20});
lastcond=subs(lastcond,{x1,x2},{x1-20,x2-20});

toc
