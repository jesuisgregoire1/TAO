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
Co = ctrb(A,B);
Obs = obsv(A,C);
n = size(A,2);
m = size(B,2);

if rank(Co) == n(1) && rank(Obs) == n(1)
    disp('Sistemul este controlabil si observabil => REALIZARE MINIMALA');
else
    disp('NU ESTE O REALIZARE MINIMALA');
end

Q = eye(n);
R = eye(m);
%%
X = icare(A, B, Q, R);
F = -R \ B' * X;

disp('Valorile proprii ale sistemului in bucla inchisa:');
disp(eig(A + B*F));
disp('Valorile proprii ale solutiei Riccati:');
disp(eig(X));
% valorile proprii ale sist in bucla inchisa au partea reala < 0, asadar
% sistemul este stabilizat
% val pr ale solutiei X au partea reala > 0, deci sol este pozitiv definita

%% Metoda Vectorilor Schur
Hf = [A -B/R*B' ;
     -Q   -A'   ];
 
[U, S] = schur(Hf);
[U, S] = ordschur(U, S, 'lhp'); %se ordoneaza matricile U si S a.i. val pr stabile
                                %ale lui A sa fie asezate in partea stanga
                                %(coloanele din stanga)
X_schur = U(n+1:2*n, 1:n) / U(1:n, 1:n); %solutia
K_schur = -R \ B' * X_schur;             %reactia

disp('Norma diferentei dintre sol obtinuta cu MVS si sol obtinuta cu care():');
disp(norm(X - X_schur));
% cele doua solutii sunt foarte apropiate (diferenta e de ordinul e-15)

disp('MVS - nr de conditionare:');
disp(cond(U(1:n, 1:n)));
% numarul de cnditionare e apropiat de 1, ceea ce indica o buna
% conditionare numerica a matricei U(1:n, 1:n), deci operatia de inversare a 
% ei nu este sensibila la schimbari mici in date, prin urmare rezultatul
% este precis
