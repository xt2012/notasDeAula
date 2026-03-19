prt =false;

function TransformadaComb(N,T,cor)
    w=linspace(-3* pi  ,3* pi  ,100000);
    C=0;
    for n=-N:N     C = C+ exp(- i *w*n*T) ;  end
    plot(w,abs(C),cor);
    grid on;
endfunction

function [w,Ys]=TF_fga(f,wmax,a,b)
    ws = 40*wmax; %amostar com taxa muito alta
    Ts=2* pi  /ws;
    Ns=floor((b-a)/Ts)+1;
    t=linspace(a,b,Ns);
    ys=f(t); % superamostrar sinal continuo com ws=40*wmax
    Nw=3600 ;%Calcular 1 periodo da tdft com muitos pontos
    w=linspace(-wmax,wmax,Nw);
    w_n=2* pi  *(w/ws);
    n=[0:Ns-1];
    Ys=(b-a)/Ns*ys*exp(- i *n'*w_n);
endfunction

function PlotTF(f,wmax,a,b,cor)
   [w,Ys]=TF_fga(f,wmax,a,b);
   subplot(121)
    title("f(t) continuo e nao periodico");
    xlabel("t (s)"); ylabel("f(t)");
    t=linspace(a,b,10000);
    ys=f(t);
    plot(t,ys,cor);
 %   axis([a,b,min(ys)-0.1*abs(min(ys)),max(ys)+0.1*abs(max(ys))]);
    hold on; grid on;
   subplot(122)
    title("F(jw) nao periodica e continua");
    xlabel("w (rad/s)"); ylabel("|F(jw)|");
    plot(w,abs(Ys)/max(abs(Ys)),cor)
    hold on; grid()
endfunction

