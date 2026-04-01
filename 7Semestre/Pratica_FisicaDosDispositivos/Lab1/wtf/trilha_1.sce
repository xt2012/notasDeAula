//Programa:trilha_1.sce
//Ajuste de Modelo Linear da Trilha-1 por Mínimos Quadrados do Exp#1
clear;

//Pontos Experimentais
N = 16;
//Coloque aqui os 16 Valores Experimentais medidos na Trilha-1
Vi = [ 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ]; //[V]
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
xstring(0.1,0.8,['EQM = ' string(EQM)]);


