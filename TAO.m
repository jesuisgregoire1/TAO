clc;
clear;
%Transfer Function
% ia = 0.2
n1 = 102.08;% n1 = nof teeth gear1 = 102.08
n2 = 1;% n2 = nof teeth gear2 = 1
T = 42 * 9.81 * 10^(-3);% T = torque developed by the motor Nm = 42 * 9.81
ia = 0.2;
K = T / ia;% K = mototr torque constant = T / ia
J1 = 3.5 * 10^-6;% J1 = moment of inertia of the rotor of the motor kg/m2
% = 3.5 x 10^-6 kg*m^2
%J2 = 0.15625;
J2 = 0.15625 * 10^-3;
% J2 = moment of inertia of the load kg/m2
% = 0.15625kg*m^2
ea = 12;% ea = applied armature voltage, V = 12
w = 10.472;% angular velocity = 10.472 rad/s
Kb = 0.014;% Kb = 0.014 V/(rad/s)
eb = Kb * w;% eb =  back emf, V
Ra = (ea - eb)/ia;
%%
% Transfer Function
num = (n1/n2) * K;
den = [Ra*(J1 + (n1/n2)^2 * J2), K*Kb, 0];
transfer_function = tf(num, den);
[A, B, C, D] = tf2ss(num, den);
Co = ctrb(A, B);
Obs = obsv(A, C);
n = size(A, 2);
m = size(B, 2);

if rank(Co) == n(1) && rank(Obs) == n(1)
    disp('Sistemul este controlabil si observabil => REALIZARE MINIMALA');
else
    disp('NU ESTE O REALIZARE MINIMALA');
end

Q = eye(n) * 1* 0.1;
R = eye(m) * 1 * 100000;
t = 0:1:1000;  
sys = ss(A, B, C, D);
input_signal = 1 * ones(size(t));
[ys, t, xs] = lsim(sys, input_signal, t);
figure;
plot(t, ys);
title('Comportamentul sistemului caz nominal de functionare');
xlabel('Timp');
ylabel('Iesirea sistemului');
grid on;

%%
X = icare(A, B, Q, R);
F = -R \ B' * X;

disp('Valorile proprii ale sistemului in bucla inchisa:');
disp(eig(A + B*F));
%valorile proprii au partea reala < 0 => sistemul este stabil
disp('Valorile proprii ale solutiei Riccati:');
disp(eig(X));
% val pr ale solutiei X au partea reala > 0, deci sol este pozitiv definita

%% Metoda Vectorilor Schur
% Definirea matricei Hf
Hf = [A -B / R * B';
      -Q   -A'];

% Decompunerea Schur
[U, S] = schur(Hf);
[U, S] = ordschur(U, S, 'lhp'); % Ordonarea matricilor U și S astfel încât valorile proprii stabile ale lui A să fie în partea stângă

% Soluția
X_schur = U(n + 1:2 * n, 1:n) / U(1:n, 1:n); 
K_schur = R \ B' * X_schur;

% Calculul normei diferenței între soluțiile obținute cu metoda Schur și `care()`
disp('Norma diferenței între soluțiile obținute cu metoda Schur și `care()`:');
%disp(norm(X - X_schur));

% Numărul de condiționare al lui U
%disp('Numărul de condiționare al lui U');
%disp(cond(U(1:n, 1:n))); % < 10 => Bine condiționată

% Verificarea stabilității sistemului în buclă închisă
disp('Valorile proprii ale sistemului în buclă închisă < 0 => STABIL');
disp(eig(A - B * K_schur)); 

% Plotarea răspunsului inițial al sistemului în buclă închisă
[yn, t, xn] = initial(ss(A - B * K_schur, zeros(size(B)), [eye(n); -K_schur], zeros(n + size(K_schur, 1), size(B, 2))), [1; 0]);
% Subgraficul pentru prima coloană a ieșirii
subplot(3, 1, 1);
plot(t, yn(:, 1), 'LineWidth', 2);
title('Răspunsul Inițial al Sistemului în Buclă Închisă');
xlabel('Timp');
ylabel('Acceleratia unghiulara');
grid on;
% Subgraficul pentru a doua coloană a ieșirii
subplot(3, 1, 2);
plot(t, yn(:, 2), 'LineWidth', 2);
%title('');
xlabel('Timp');
ylabel('Viteza unghiulara');
grid on;
% Subgraficul pentru a treia coloană a ieșirii
subplot(3, 1, 3);
plot(t, yn(:, 3), 'LineWidth', 2);
%title('Răspunsul Inițial al Sistemului în Buclă Închisă');
xlabel('Timp');
ylabel('Comanda');
grid on;


%% Metoda Vectorilor Schur in caz de avarie

Ra_avarie = Ra * 10;  % Dublarea rezistenței
% Recalculare parametrii
num_avarie = (n1/n2) * K;
%den_avarie = [Ra_avarie*(J1 + (n1/n2)^2 * J2), K*Kb, 0];
%den_avarie = [100, -1.5, 0]; 
den_avarie = [-100, -1.5,0];

% Creare sistem în spațiul stărilor pentru cazul avariat
[Af, Bf, Cf, Df] = tf2ss(num_avarie, den_avarie);
sysf = ss(Af, Bf, Cf, Df);


%% Metoda Vectorilor Schur
Hf = [Af -Bf/R*Bf' ;
     -Q   -Af'   ];
 
[U, S] = schur(Hf);
[U, S] = ordschur(U, S, 'lhp'); 
                                
X_schurf = U(n+1:2*n, 1:n) / U(1:n, 1:n); 
K_schurf = R \ Bf' * X_schurf;            

disp('MVS - nr de conditionare:');
disp(cond(U(1:n, 1:n))); % < 10 => Bine conditionata

disp('Valorile priprii ale sistemului in bucla inchisa < 0 => STABIL')
disp(eig(Af-Bf*K_schurf));%valorile proprii au partea reala < 0 => sistemul este stabil
[yf, t, xf] = initial(ss(Af-Bf*K_schurf, zeros(size(Bf)), [eye(n); -K_schurf], zeros(n+size(K_schurf,1), size(Bf,2))), [1;0]);

figure;
% Subgraficul pentru prima coloană a ieșirii
subplot(3, 1, 1);
plot(t, yn(:, 1), 'LineWidth', 1, 'DisplayName', 'Nominal');
hold on;
plot(t, yf(:, 1), 'LineWidth', 1, 'DisplayName', 'Avariat');
title('Acceleratia unghiulara');
xlabel('Timp');
ylabel('Valoare');
legend;
grid on;

% Subgraficul pentru a doua coloană a ieșirii
subplot(3, 1, 2);
plot(t, yn(:, 2), 'LineWidth', 1, 'DisplayName', 'Nominal');
hold on;
plot(t, yf(:, 2), 'LineWidth', 1, 'DisplayName', 'Avariat');
title('Viteza unghiulara');
xlabel('Timp');
ylabel('Valoare');
legend;
grid on;

% Subgraficul pentru a treia coloană a ieșirii
subplot(3, 1, 3);
plot(t, yn(:, 3), 'LineWidth', 1, 'DisplayName', 'Nominal');
hold on;
plot(t, yf(:, 3), 'LineWidth', 1, 'DisplayName', 'Avariat');
title('Comanda');
xlabel('Timp');
ylabel('Valoare');
legend;
grid on;

