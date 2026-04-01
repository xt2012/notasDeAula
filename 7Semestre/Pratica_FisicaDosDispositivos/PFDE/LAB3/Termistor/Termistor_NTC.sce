//%Programa: Termistor_NTC.sce
//Cálculo da temperatura do termistor NTC 5D-9
clear;
R2 = 10; //Coloque aqui o valor do Resistor R2 [ohms]
Tamb = 21; //Temperatura Ambiente [Graus Centígrados]
K = 11e-3; //[W/Graus] Coeficiente Térmico de Dissipação por Convecção
//Pontos Experimentais
N = 8;
//Coloque aqui os valores exatos medidos de Va e Vb:
Va = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]; //[V]
Vb = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]; //[V]
//Cáculo da Corrente (I)) e da Resistência (R1)
I = Vb./R2; //[A]
Vab = Va - Vb; //[V]
Pe = Vab .* I; // [W] Potência Elétrica sobre o termistor
R1 = Vab./I; //[ohm] Resistência do Termistor
//Sabendo que:Pe = P = K*(T-Tamb); //[W] Potência Térmica Dissipada
T = (Pe./K)+Tamb; //[Graus] Temperatura do Termistor
scf(1);
plot(T, R1,'ro');
plot(T, R1,'b');
xlabel('T[Graus]');
ylabel('R1[Ohms]');
title('Gráfico: R1 vs. T');

