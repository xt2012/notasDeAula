//Exemplo: Simula Resposta ao Degrau do Filamento
//%Programa:Simula_Degrau.sce
clear;
//Parâmetros
M = 4.6e-3; //[g] Massa do Filamento
R2 = 1; //[ohm]
Klt = 0.011; //[ohm/K] Constante do Filamento (Resistência)
Kr = 3.3e-14; //[W/K^4] Constante do Filamento (Radiação)
Tamb = 300;//[K] Temperatura Inicial (Ambiente)
cp = 0.13;//[J/g.K] Calor Específico
Ct = cp*M; //[J/k] Capacidade de Calor 
tmax = 800e-3; //[s] Tempo da Simulação
N = 10000; //Número de Pontos da Simulação
t = linspace(0,tmax,N); //Base de Plotagem do Mod
T = 0*t; //[K] Temperatura do Filamento
R1 = 0*t; //[ohms] Resistência da Lâmpada
I = 0*t; //[A] Corrente
dt = t(2)-t(1); //Passo de Tempo [s]
V = 12;//[V] Tensão Aplicada
k = 1; //Índice do Passo Inicial
T(1) = 300; //Temperatura Inicial
R1(1)= 3.4; //[ohms]
I(1)= V/(R1(1)+R2); //[A]
//LOOP DE SIMULAÇÃO
while(k<N)
Pe = Klt*T(k)*((V/(Klt*T(k)+R2))^2); // [watts] Potência Elétrica
Pr = Kr*((T(k))^4); //[watts] Potência Irradiada (Corpo Negro)
T(k+1) = T(k) + (dt/Ct)*(Pe - Pr); //Cálculo da Evolução da Temperatura
R1(k+1) = Klt*T(k+1); //[ohms] Resistência da Lâmpada
I(k+1)= V/(R1(k+1)+R2); //[A] Corrente na Lâmpada
k = k + 1; //Incrementa o Passo do tempo discreto   
end
//Estimando o Tempo de Descida (Critério de I1:10% e I2:90%)
Imax = max(I); //[A]
Imin = min(I);
dI = Imax - Imin;
I1 = I(min(find(I<(Imin + 0.9*dI)))); //[K]
t1 = t(min(find(I<(Imin + 0.9*dI)))); //[s]
scf(2);
plot(0,max(I),'bo');
plot(max(t),min(I),'bo');
plot(t1,I1,'ro');
plot(t1,I1,t1,0,'r');
plot([t1 t1],[I1 0],'r');
I2 = I(min(find(I<(Imin + 0.1*dI)))); //[K]
t2 = t(min(find(I<(Imin + 0.1*dI)))); //[s]
scf(2);
plot(t2,I2,'ro');
plot([t2 t2],[I2 0],'r');
//Tempo de Descida:
Td = t2 - t1; //[s]
mprintf('Tempo de Descida: Td = %.1e [s]\n',Td);
//GRÁFICOS:
scf(1);
plot(t,T,'r');
title('Modelo Dinâmico: Lâmpada Elétrica de Filamento');
xlabel('t [s]');
ylabel('T [K]');
scf(2);
plot(t,I,'b');
title('Modelo Dinâmico: Lâmpada Elétrica de Filamento');
xlabel('t [s]');
ylabel('I [A]');
scf(3);
plot(t,R1,'g');
title('Modelo Dinâmico: Lâmpada Elétrica de Filamento');
xlabel('t [s]');
ylabel('R1 [ohm]');
