clc;
clear;
%Transfer Function
%%
% n1= nof teeth gear1
% n2= nof teeth gear2
% K = mototr torque constant
% Ra = 
% J1 = moment of inertia of the rotor of the motor kg/m2
% J2 = moment of inertia of the load kg/m2
% Kb =
% T = torque developed by the motor Nm
% ea = applied armature voltage, V
% eb =  back emf, V

%%
ia = T / K;
Ra = (ea - eb)/ia;
Kb = eb/delta_theta1;

%TransferFunction
num = (n1/n2) * K;
den = [Ra*(J1 + (n1/n2)^2 * J2), K*Kb];

transfer_function = tf(num,den); %for testing purposes
[A,B,C,D] = tf2ss(num,den);