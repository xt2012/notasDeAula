//Programa: LED_LDR_MQ_1.sce
//Ajuste do Modelo Não-linear do LED e LDR 
clear;
//Dados Experimentais - 1
N = 6;
//Descartando ponto com Iled = 0, onde o LED está apagado.
R2 = 979; //[ohms]
R1=101.3;
VR2 = [0.499 1.002 1.501 2.01 2.50 3.00]; //[V]
VR1 = [0.985 1.424 1.725 1.951 2.14 2.3]
Iled = VR2./R2;//[A]
Rldr = [1.115e3 0.748e3 0.601e3 0.518e3 0.464e3 0.424e3];//[ohms]
scf(1);
plot(Iled,Rldr,'or');
ylabel('Rldr [ohms]');
xlabel('Iled [A]');
title('Resistência do LDR vs. Corrente do LED');
//MODELO: 
//Gldr = C1*sqrt(Iled) + C2
Rldr2 = VR1./R1;
Gldr = 1./Rldr2;
//Funções Base Não-Ortogonais
g1 = sqrt(Iled);
g2 = ones(1,N);
a11 = sum(g1.*g1);
a12 = sum(g1.*g2);
a21 = sum(g2.*g1);
a22 = sum(g2.*g2);
b1 = sum(Gldr.*g1);
b2 = sum(Gldr.*g2);
A = [ a11 a12; a21 a22 ];
b = [b1; b2];
C = A\b;
disp('C = ',C);
M = 100;
Iledc = linspace(min(Iled),max(Iled),M);
Gldrc = C(1)*sqrt(Iledc) + C(2);
Rldrc = 1./Gldrc;
plot(Iledc,Rldrc,'b');
Gldr0 = C(1)*sqrt(Iled) + C(2);
Rldr0 = 1./Gldr0;
//Erro Quadrático Médio Percentual
EQMP = 100*(1/N)*sum(((Rldr - Rldr0)./Rldr).^2);
disp('EQM% = ',EQMP);






