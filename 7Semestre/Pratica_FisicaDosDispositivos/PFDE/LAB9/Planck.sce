//%Programa:Planck.sce
clear;
h = 6.62e-34; //[J.s]
e = 1.602e-19; //[C]
c = 2.99e8; //[m/s]
lambda = [ 430 565 590 627 700 850 ]; //[nm]
Vs = [2.50 1.85 1.70 1.65 1.73 1.20]; //[V]
//Modelo: lambda = K * (1/Vd)
g = 1./Vs; //função-base
f = lambda; 
K = sum(f.*g)/sum(g.*g); //Ajuste da Constante 
hexp = (1e-9)*K*e/c; //Constante de Planck [J.s]
mprintf('h(experimental) = %.2e [J.s]\n',hexp);
scf(1);
N = 100;
V = linspace(min(Vs),max(Vs),N);
Li = (1e9)*(h*c/e)*(1./V); //Ideal [nm]
La = K * (1./V); //Ajustado [nm]
plot(V, Li,'b'); //Ideal (azul)
plot(V, La,'r'); //Ajustado (vermelho)
plot(Vs, lambda,'or'); //Experimental
title('Lambda vs. Vs');
xlabel('Vs [V]');
ylabel('Lambda [nm]');










