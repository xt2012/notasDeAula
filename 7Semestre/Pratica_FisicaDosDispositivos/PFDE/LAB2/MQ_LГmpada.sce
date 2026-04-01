//Exemplo: Ajuste de Modelo da Lâmpada por Mínimos Quadrados
//%Programa:MQ_Lâmpada.sce
clear;

//Pontos Experimentais
N = 7; //Número de pontos experimentais
//SUBSTITUA AQUI OS VALORES MEDIDOS:
yp = [0, 0.12, 0.15, 0.19, 0.21, 0.23, 0.26]; //Corrente [I]
xp = [0, 1.87, 3.81, 5.75, 7.7, 9.66, 11.62]; //Tensão [V]


plot(xp,yp,'or');
R1 = xp./yp;
//Modelo Adotaldo: y = k*(x^p)
p = 3/5; //Expoente Fracioqnário (Modelo)
g = xp.^p; //Função Base
k = sum(yp.*g)/sum(g.*g); //Ajuste da Constante do Modelo

xc = linspace(min(xp),max(xp),100); //Base de Plotagem do Modelo
yc = k*(xc.^p); //Modelo Ajustado
//plot(xc,yc,'b');
title('Modelo Ajustado (azul) e Pontos Experimentais (vermelho)');

//Erro Quadrático Médio
ym = k*(xp.^p)//Valores da corrente a partir do modelo ajustado
EQM = (1/N)*sum((ym-yp).^2)
disp(EQM);
Q2 = abs(k*((-5)^p));
disp(Q2)

