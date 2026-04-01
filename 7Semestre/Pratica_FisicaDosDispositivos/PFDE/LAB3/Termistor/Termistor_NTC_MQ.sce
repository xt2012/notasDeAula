//Programa: Termistor_NTC_MQ.sce
//Ajuste por Mínimos Quadrados
//Modelo de Steinhart-Hart
clear;

//Dados do Fabricante para o Termistor NTC 10D-9
R0 = 5; //[ohms] //Resistência Aproximada na Temperatura Ambiente
K = 11e-3; //[W/K] Coeficiente da Convecção Natural (Newton)
T0 = 23; //[Graus] Temperatura Ambiente Medida
Imax = 2; //[A] Corrente Máxima
Rmax = 0.458; //[Ohms] Resistência na Corrente Máxima
Pmax = Rmax*(Imax^2); //[W] Potência Máxima
Tmax = Pmax/K + T0; //[Graus] Temperatura Máxima
//Cálculo do Coeficiente B do Modelo de Steinhart-Hart
B = (1/(1/Tmax - 1/T0))*log(Rmax/R0); // [Graus]
//Cálculo de rinf do Modelo de Steinhart-Hart  
rinf = R0*exp(-B/T0); //[Ohms]

//Resistor
R2 = 9.6; //[ohms]

//Pontos Experimentais
N = 8;
Va = [0.503 0.987 1.506 2.0 2.5 3.0 3.49 4.01]; //[V]
Vb = [0.311 0.614 0.976 1.374 1.807 2.26 2.77 3.29]; //[V]

//Cáculo da Corrente (I) e da Resistência (R1)
I = Vb./R2; //[A] Corrente no Termistor
Vab = Va - Vb; //[V] Tensão no Termistor
Pe = Vab .* I; //[W] Potência Elétrica Dissipada
R1 = Vab./I -1; //[ohm] Resistência do Termistor
P = Pe; //[W] Equilíbrio Térmico
T = (P/K) + T0; //[K] Convecção Natural (Newton)

scf(1);
plot(T,R1,'or'); //Pontos vermelhos experimentais
xlabel('T [Graus]');
ylabel('R1 [Ohms]');
title('Gráfico: R1 vs. T');

//AJUSTE DO MODELO DE STEIHART-HART
//R1 = rinf * exp(B/T)
//Pseudo-linearização:
//log(R1) = log(rinf) + log(exp(B/T)) = C + B/T
//Onde: C = alfa(1) e B = alfa(2)

//Funções Base Não-Ortogonais:
g1 = ones(1,N);
g2 = 1./T;

//Coeficientes da Matriz A:
a11 = sum(g1.*g1);
a12 = sum(g1.*g2);
a21 = sum(g2.*g1);
a22 = sum(g2.*g2);

//Coeficientes do Vetor b:
b1 = sum(log(R1).*g1);
b2 = sum(log(R1).*g2);

//Sistema:
A = [ a11 a12; a21 a22 ];
b = [b1; b2];

//Solução do Sistema de Equações
alfa = A\b;

//Plota curva do modelo ajustado
M = 100;
Tc = linspace(min(T),max(T),M);
gc1 = ones(1,M);
gc2 = 1./Tc;
R1c = exp(alfa(1))*exp(alfa(2)*gc2);
plot(Tc,R1c,'b'); //Curva azul

//Erro Quadrático Médio (EQM)
gc1 = ones(1,N);
gc2 = 1./T;
R1p = exp(alfa(1))*exp(alfa(2)*gc2);
EQM = (1/N)*sum((R1 - R1p).^2);
mprintf('EQM = %.1e \n',EQM);

//Curva obtida com dados do fabricante
B = T0*log(R0/rinf);
R1f = rinf.*exp(B./Tc);
plot(Tc,R1f,'r'); //Curva vermelha





