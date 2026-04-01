//Programa: Lab_10.sce
clear;
//Pontos Experimentais (Item-1, +10V)
Vin = 2; //[Vp-p]
Vout = [36e-3 290e-3 1.74 2.00 2.00]; //[Vp-p]
Vout_inv = [2.72e-3 2.16e-3 2.00e-3 3.6e-3 15.6e-3]; //[Vp-p]
f = [100 1e3 1e4 1e5 1e6];//[Hz]
w = 2*%pi*f; //[rad/s]
Hp = Vout/Vin;
Hp_inv = Vout_inv/Vin;
atten = 20*log(Hp)
atten_inv=20*log(Hp_inv)

//Modelo:
C = 0.1e-6;//[F]
R = 470;//[ohms]
Reff = R/2;//[ohms]
x = abs((%i*w*C*Reff)./(1+%i*w*C*Reff)); //Função Base
y = Hp; 
k = sum(y.*x)/sum(x.*x); //Ajuste da Constante do Modelo
M = 100; //Número de pontos para plotar modelo
fp = logspace(log10(100),log10(1e6),M); //[Hz]
wp = 2*%pi*fp; //[rad/s]
//Função:
Hm = k*abs((%i*wp*C*Reff)./(1+%i*wp*C*Reff));

//Gráfico: |H(jw)| Log-Linear
scf(1);
plot2d("ln",wp,Hm);
plot(w,Hp,'or');
xlabel('w[rad/s]');
ylabel('|H(jw)|');
title('Gráfico: |H(jw)| vs. w');
