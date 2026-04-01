////%Programa:Optoacoplador.sce
//clear;
//R1 = 0.996e3; //[ohms]
//R2 = 0.99e6; //[ohms]
//
////Parte-I: Fotocondutivo
//Vdc1a = [0 5 10 15 20 25];
//Vdc2a = 12.0; //[V]
//V1a = [0 3.18 8.04 12.96 17.85 22.5];
//V2a = [0 0.422 1.076 1.632 2.08 2.37];
//Vled1a = Vdc1a - V1a;
//Vled2a = Vdc2a - V2a;
//Iled1a = V1a./R1;
//Iled2a = V2a./R2;
//P1a = Vled1a.*Iled1a;
//P2a = Vled2a.*Iled2a;
//scf(1);
//plot(P1a,P2a,'or');
//
////Parte-II: Fotovoltaico
//Vdc1b = [0 5.0 10.00 15.00 20.0 25];
////Vdc2b = 12.05; //[V]
//V1b = [0 3.2 8.0 12.95 17.8 22.9];
//V2b = [0 -0.372 -0.902 -1.291 -1.404 -1.435];
//Vled1b = Vdc1b - V1b;
////Vled2b = Vdc2b - V2b;
//Vled2b = V2b;
//Iled1b = V1b./R1;
//Iled2b = V2b./R2;
//P1b = Vled1b.*Iled1b;
//P2b = Vled2b.*Iled2b;
//scf(2);
//plot(P1b,P2b,'or');

//Parte-III: I vs V Fotovoltaico
// Com Rmult = 20e6 => Vca = 1258e-3
R = [19.99e3 0.498e6 1.001e6 1.499e6 2e6 2.5e6 3e6 3.5e6 4.01e6 4.51e6 5e6 9.96e6];
V = [-38.8e-3 -885e-3 -1.404 -1.460 -1.473 -1.477 -1.490 -1.493 -1.494 -1.494 -1.499 -1.504];
Rmult = 20e6;//20e6;
Reff = 1./((1./R)+(1./Rmult));
Reff(12) = 20e6; //Correção: Aqui apenas R do multímetro.
I = V./Reff;
IL = I(1); //IL ~ I(R pequeno), V ~ 0
scf(3);
plot(V,I,'ro');
plot(V,I,'b');
//replot([0 0 1.4 6e-7],'tight');
title('I vs. V (Modo Fotovoltaico)');
xlabel('V [volt]');
ylabel('I [A]');
//
////Ajuste Linear da Parte-I
////Modelo: y = eta*x 
//g = P1a; //Função Base
//f = P2a; 
//eta = sum(f.*g)/sum(g.*g); //Ajuste da Constante do Modelo
//disp(eta);
//scf(1);
//N = 100;
//x = linspace(0,max(P1a),N);
//y = eta*x;
//plot(x, y,'b');
//title('P2 vs. P1 (Fotocondutivo)');
//xlabel('P1 [W]');
//ylabel('P2 [W]');

////Ajuste Linear da Parte-II
////Modelo: y = eta*x 
//g = P1b; //Função Base
//f = P2b; 
//eta = sum(f.*g)/sum(g.*g); //Ajuste da Constante do Modelo
//disp(eta);
//scf(2);
//N = 100;
//x = linspace(0,max(P1b),N);
//y = eta*x;
//plot(x, y,'b');
//title('P2 vs. P1 (Fotovoltaico)');
//xlabel('P1 [W]');
//ylabel('P2 [W]');










