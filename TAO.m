clc;
clear;
%Transfer Function
%%
% n1 = nof teeth gear1 = 102.08
% n2 = nof teeth gear2 = 1
% K = mototr torque constant = T / ia
% Ra = 
% J1 = moment of inertia of the rotor of the motor kg/m2
% = 3.5 x 10^-6 kg*m^2
% J2 = moment of inertia of the load kg/m2
% = 0.15625kg*m^2
% Kb = 0.014 V/(rad/s)
% T = torque developed by the motor Nm = 42 * 9.81
% ea = applied armature voltage, V = 12
% eb =  back emf, V
% angular velocity = 10.472 rad/s
% ia = 0.2
%% 
n1 = 102.08;
n2 = 1;
T = 42 * 9.81 * 10^(-3);
ia = 0.2;
K = T / ia;
J1 = 3.5 * 10^-6;
%J2 = 0.15625;
J2 = 0.15625 * 10^-3;
ea = 12;
w = 10.472;
Kb = 0.014;
eb = Kb * w;
Ra = (ea - eb)/ia;
%%

%TransferFunction
num = (n1/n2) * K;
den = [Ra*(J1 + (n1/n2)^2 * J2), K*Kb, 0];

transfer_function = tf(num,den); %for testing purposes
[A,B,C,D] = tf2ss(num,den);