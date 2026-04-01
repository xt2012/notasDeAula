//%Programa: Termistor_NTC.sce
//Cálculo da temperatura do termistor NTC 5D-9
clear;
close;
R2 = 9.6; //Coloque aqui o valor do Resistor R2 [ohms]
Tamb = 21; //Temperatura Ambiente [Graus Centígrados]
K = 11e-3; //[W/Graus] Coeficiente Térmico de Dissipação por Convecção
//Pontos Experimentais
N = 8;
//Coloque aqui os valores exatos medidos de Va e Vb:
Va = [0.503 0.987 1.506 2.0 2.5 3.0 3.49 4.01]; //[V]
Vb = [0.311 0.614 0.976 1.374 1.807 2.26 2.77 3.29]; //[V]
//Cáculo da Corrente (I)) e da Resistência (R1)
I = Vb./R2; //[A]
Vab = Va - Vb; //[V]
Pe = Vab .* I; // [W] Potência Elétrica sobre o termistor
R1 = Vab./I -1; //[ohm] Resistência do Termistor
//Sabendo que:Pe = P = K*(T-Tamb); //[W] Potência Térmica Dissipada
T = (Pe./K)+Tamb; //[Graus] Temperatura do Termistor
scf(1);
plot(T,R1,'ro');
plot(T,R1,'b');
xlabel('T[Graus]');
ylabel('R1[Ohms]');
title('Gráfico: R1 vs. T');

//Modelo Adotaldo: y = k*(x^p)
r00=0.00018;
R0=5;
T0=298.15;
B = T0*log10(R0/r00);
Tcal = (B ./ (log10(R1/r00))) - T0 + Tamb;

scf(2);
plot(Tcal, R1, 'r');
plot(T, R1, 'bo');
xlabel('T[graus]');
ylabel('R1[ohms]');
title('Gráfico: R1 vs T Teórico');
