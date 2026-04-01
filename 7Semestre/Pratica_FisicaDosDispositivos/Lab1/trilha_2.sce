//Programa:trilha_1.sce
//Ajuste de Modelo Linear da Trilha-1 por Mínimos Quadrados do Exp#1
clear;

//Pontos Experimentais
N = 8;
//Coloque aqui os 16 Valores Experimentais medidos na Trilha-1
Vi = [ 0.5 0.75 1.46 2.01 2.68 3.34 4.05 4.77 ]; //[V]
i = 1:N;
W = 0.01; //[m]
xi = W.*(i-1)+W/2;
plot(xi,Vi,'or');

//Modelo Linear Adotado: y = k*(x - W/2);
g = xi - W/2; //Função Base
k = sum(Vi.*g)/sum(g.*g); //Ajuste da Constante do Modelo

xc = linspace(min(xi),max(xi),100);//Base de Plotagem
yc = k*(xc - W/2); //Modelo Ajustado
plot(xc,yc,'b');
title('Trilha-1');
ylabel('Potencial [V]');
xlabel('Posição [m]');

//Erro Quadrático Médio
EQM = (1/N)*sum((Vi-k*(xi-W/2)).^2);
xstring(0.02,4,['EQM = ' string(EQM)]);


