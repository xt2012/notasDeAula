//Programa:trilha_1.sce
//Ajuste de Modelo Linear da Trilha-1 por Mínimos Quadrados do Exp#1
clear;

//Pontos Experimentais
N = 16;
//Coloque aqui os 16 Valores Experimentais medidos na Trilha-1
Vi = [ 0 0.201 0.314 0.429 0.544 0.641 0.77 0.899 1.04 1.17 1.273 1.385 1.501 1.668 1.885 2.03 ]; //[V]
//Vi = [ 0 0.405 0.711 1.058 1.396 1.618 1.805 1.97 ]; //[V]
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


