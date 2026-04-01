//%Programa:Termoeletrico.sce
clear;
//Parte-I: V2 vs. P1
V1a = [0 1.00 2.0 3.0 4.0 5.0 6.0 7.0];
I1a = [0 0.22 0.44 0.63 0.85 1.02 1.15 1.32];
V2a = [0 3.8e-3 5.1e-3 8.7e-3 9.7e-3 12.1e-3 13.6e-3 15.4e-3];
P1a = V1a.*I1a;
scf(1);
plot(P1a, V2a,'or');

//Parte-II: P2 vs. P1
V1b = [0 1.00 2.0 3.0 4.0 5.0 6.0 7.0];
I1b = [0 0.2 0.41 0.62 0.82 0.95 1.07 1.21];
V2b = [0 1.4e-3 2.7e-3 3.6e-3 6.7e-3 10.0e-3 10.7e-3 13.0e-3];
R2b = 2.2e3; //[ohms]
I2b = V2b./R2b;
P1b = V1b.*I1b;
P2b = V2b.*I2b;
scf(2);
plot(P1b,P2b,'ob');

//Ajuste Linearizado da Parte-I
//Modelo: y = k*sqrt(x) 
g = sqrt(P1a); //Função Base
f = V2a; 
k = sum(f.*g)/sum(g.*g); //Ajuste da Constante do Modelo
disp(k);
scf(1);
N = 100;
x = linspace(0,max(P1a),N);
y = k*sqrt(x);
plot(x, y,'b');
title('V2 vs. P1');
xlabel('P1 [W]');
ylabel('V2 [V]');

//Ajuste Linear da Parte-II
//Modelo: y = eta*x 
g = P1b; //Função Base
f = P2b; 
eta = sum(f.*g)/sum(g.*g); //Ajuste da Constante do Modelo
disp(eta);
scf(2);
N = 100;
x = linspace(0,max(P1b),N);
y = eta*x;
plot(x, y,'r');
title('P2 vs. P1');
xlabel('P1 [W]');
ylabel('P2 [W]');









