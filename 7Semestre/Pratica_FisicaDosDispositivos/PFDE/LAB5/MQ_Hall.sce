//%Programa:MQ_Hall.sce
//Exemplo: Ajuste de Modelo por Mínimos Quadrados do Sensor Hall
clear;
//Dados Experimentais:
Vs_max = 4.89; //[V] Valor medido no ângulo teta = 0 graus.
ks = 3.125e-3; //[V/Gauss] Constante do Sensor (Dados do Fabricante)
N = 12; //Número de Pontos Experimentais
Vso = 2.51; //[V] Valor medido na ausência de campo magnético.
teta = [0, 30,60,90,120,150,180,210,240,270,300,330]; //ângulo [graus]
//Valores medidos da tensão de saída Vs para cada valor de`ângulo.
//COLOQUE AQUI AS SUAS MEDIDAS!
Vs = [4.43, 3.78, 2.55, 1.657, 1.228, 0.860, 1.070, 1.618, 2    .25, 3.18, 4.48]
//Cálculo de B0:
B0 = (Vs_max - Vso)/ks; //[Gauss]
//Cálculo de Bz em função de B0, para cada ângulo:
Bz = B0*cos((2*%pi/360)*teta); //[Gauss]

//Ajuste do Modelo Não-linear
g1 = cos((2*%pi/360)*teta);
g2 = ones(1,N);
a11 = sum(g1.*g1);
a12 = sum(g1.*g2);
a21 = sum(g2.*g1);
a22 = sum(g2.*g2);
b1 = sum(Vs.*g1);
b2 = sum(Vs.*g2);
A = [ a11 a12; a21 a22];
b = [b1; b2];
alfa = A\b;
disp('alfa (Não-linear) = ',alfa);
M = 100;
angulo = linspace(min(teta),max(teta),M);
gc1 = cos((2*%pi/360)*angulo);
gc2 = ones(1,M);
Vs_lin = alfa(1)*gc1 + alfa(2)*gc2; 

//Plotando o Modelo Não-linear
scf(1);
plot(teta,Vs,'or'); //Pontos Experimentais
plot(angulo,Vs_lin,'b'); //Modelo Ajustado Não-linear
xlabel('Ângulo [graus]');
ylabel('Vs[volts]');
title('Sensor Hall');
Vs_mod_lin = alfa(1)*g1 + alfa(2)*g2;
EQM = (1/N)*sum((Vs_mod_lin - Vs).^2);
disp('EQM (Mod. Não-linear) = ',EQM); //Erro Quadrático Médio (Não-linear)

//Ajuste do Modelo Linear
g1 = Bz;
g2 = ones(1,N);
a11 = sum(g1.*g1);
a12 = sum(g1.*g2);
a21 = sum(g2.*g1);
a22 = sum(g2.*g2);
b1 = sum(Vs.*g1);
b2 = sum(Vs.*g2);
A = [ a11 a12; a21 a22];
b = [b1; b2];
alfa = A\b;
disp('alfa (Linear) = ',alfa);
M = 100;
campo_z = linspace(min(Bz),max(Bz),M);
gc1 = campo_z;
gc2 = ones(1,M);
Vs_nlin = alfa(1)*gc1 + alfa(2)*gc2; 

//Modelo Linear
scf(2);
plot(Bz,Vs,'or');
plot(campo_z,Vs_nlin,'b');
xlabel('Bz [Gauss]');
ylabel('Vs [volts]');
title('Sensor Hall');
Vs_mod_nlin = alfa(1)*g1 + alfa(2)*g2;
EQM = (1/N)*sum((Vs_mod_nlin - Vs).^2);
disp('EQM (Mod. Linear) = ',EQM); //Erro Quadrático Médio (Linear)

