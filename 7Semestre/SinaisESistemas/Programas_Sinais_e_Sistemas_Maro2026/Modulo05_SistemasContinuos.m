prt = false;

function exemplo_A1() % resposta ao degrau
   ct=1/6.25;
   t=linspace(-ct,5*ct,1000);
   yo = @(t) (1-exp(-6.25*t)).*f_u(t);
   hold on; grid on;
   plot(t,f_u(t),'r')
   plot(t,yo(t),'b')
   axis([-ct,5*ct])
endfunction

function exemplo_A2()   % impulso
    ct=1/6.25;
    t=linspace(-ct,5*ct,3000);
    h = @(t)  6.25*exp(-6.25*t).*f_u(t);
    hold on; grid on;
    plot(t,h(t),"b")
    quiver(0,0,0,4, 'MaxHeadSize', 0.01, 'Color',
           'red', 'LineWidth', 1);
    axis([-ct,5*ct])
endfunction

function exemplo_A3o()   %seno Estado Nulo
    ct=1/6.25;
    t=linspace(-ct,20*ct,3000);
    xi = @(t)  10.5*exp(-t).*f_u(t);
    yo = @(t) (12.5*exp(-t)-12.5*exp(-6.25*t)).*f_u(t);
    hold on; grid on;
    plot(t,yo(t),"b");
    plot(t,xi(t),"r");
    axis([-ct,20*ct])
endfunction

function exemplo_A4o()   %seno Estado Nulo
   ct=1/6.25;
   t=linspace(-ct,20*ct,3000);
   xi = @(t) 10*cos(7*t).*f_u(t);
   yo = @(t) (-4.44*exp(-6.25*t)+6.66*cos(7*t-0.842)).*f_u(t);
    hold on; grid on;
    plot(t,yo(t),"b");
    plot(t,xi(t),"r");
    axis([-ct,20*ct]);
endfunction

function exemplo_A4n() %seno Entrada Nula
   ct=1/6.25;
   t=linspace(-ct,20*ct,3000);
   yn = @(t) 5*exp(-6.25*t).*f_u(t);
   hold on; grid on;
   plot(t,yn(t),"b");
   plot(t,t*0,"r");
   axis([-ct,20*ct,-10,10]);
endfunction

function exemplo_A4c() %seno Completa
   ct=1/6.25;
   t=linspace(-ct,20*ct,3000);
   xi = @(t) 10*cos(7*t).*f_u(t);
   yn = @(t) 5*exp(-6.25*t).*f_u(t);
   yo = @(t) (-4.44*exp(-6.25*t)+6.66*cos(7*t-0.842)).*f_u(t);
   yc = @(t) yn(t)+yo(t);
   hold on; grid on;
   plot(t,xi(t),"r");
   plot(t,yc(t),"b");
   axis([-ct,20*ct]);
endfunction



function y_out=Sistema_Eletrico_1aOrdem(x_in,y_ini)
    xi=[0 x_in];
    yi=[y_ini];
    for (i=2:length(xi))
          yi(i)=(1/6)*xi(i)+(5/6)*yi(i-1);
    end
    y_out=yi(2:length(xi));
endfunction

function y_out=Sistema_Eletrico_1aOrdem_T(x_in,y_ini,T)
    den=(1+6.25*T);
    a1 = 1/den;
    b0 = (6.25*T)/den;
    xi=[0 x_in];
    yi=[y_ini];
    for (i=2:length(xi))
          yi(i)=b0*xi(i)+a1*yi(i-1);
    end
    y_out=yi(2:length(xi));
endfunction


function exemplo_A5(T)
    ct=1/6.25;
    ini=0;
    fim=5*ct;
    N=ceil((fim-ini)/T +1);
    n=[-1:N-2];
    subplot(521)  %degrau
        x=[0 ones(1,N-1)];
        bar(n*T,x,0.5,'r','edgecolor','none')
        grid on;
    subplot(522) %resposta ao degrau
        y=Sistema_Eletrico_1aOrdem(x,0,T);
        bar(n*T,y,0.5,'b','edgecolor','none')
        grid on;
    subplot(523) %impulso
        x=[0 1 zeros(1,N-2)];
        bar(n*T,x,0.5,'r','edgecolor','none')
        grid on;
    subplot(524) %resposta ao impulso
        y=Sistema_Eletrico_1aOrdem(x,0,T);
        bar(n*T,y,0.5,'b','edgecolor','none')
        grid on;
    subplot(525) %exponencial
        fim=20*ct;
        N=ceil((fim-ini)/T +1);
        n=[0:N-1];
        xi = @(t) 10.5*exp(-t).*f_u(t);
        bar(n*T,xi(n*T),0.5,'r','edgecolor','none')
        grid on;
    subplot(526) % resposta a exponencial
        y=Sistema_Eletrico_1aOrdem(xi(n*T),0,T);
        bar(n*T,y,0.5,'b','edgecolor','none')
    subplot(527) %seno
        xi = @(t) 10*cos(7*t).*f_u(t);
        bar(n*T,xi(n*T),0.5,'r','edgecolor','none')
        grid on;
    subplot(528) % resposta ao seno
        y=Sistema_Eletrico_1aOrdem(xi(n*T),0,T);
        bar(n*T,y,0.5,'b','edgecolor','none')
        grid on;
    subplot(529) %seno
        bar(n*T,xi(n*T),0.5,'r','edgecolor','none')
        grid on;
    subplot(5,2,10) % resposta ao seno com cond. inicial
        y=Sistema_Eletrico_1aOrdem(xi(n*T),5,T);
        bar(n*T,y,0.5,'b','edgecolor','none')
        grid on;
endfunction

function exemplo_B1() %degrau
        ini=-0.2;
        fim=1.0;
        t=linspace(ini,fim,3000);
        hold on; grid on;
        yo= @(t) (1-1.5*exp(-30*t)+0.5*exp(-90*t)).*f_u(t);
        plot(t,f_u(t),"r--");
        plot(t,yo(t),"b");
        legend("entrada","saida");
        legend('location','southeast');
        axis([ini,fim,-0.01,1.1]);
endfunction

function exemplo_B2()   % impulso
        ini=-0.2;
        fim=1.0;
        t=linspace(ini,fim,3000);
        h= @(t) (45*exp(-30*t)-45*exp(-90*t)).*f_u(t);
		    hold on; grid on;
	 	    plot(0,0,'r')
        plot(t,h(t),"b")
		    legend("entrada","saida");
        legend('location','southeast');
        plot_arrow_fga([0,0],[0,14],'r')
endfunction

function exemplo_B3o()   %exponencial Estado Nulo
        ini=-0.2;
        fim=1.0;
        t=linspace(ini,fim,3000);
        xi= @(t) 10*exp(-15*t).*f_u(t);
        yo=@(t)(6*exp(-90*t)-30*exp(-30*t)+24*exp(-15*t) ).*f_u(t);
        hold on; grid on;
        plot(t,xi(t),"r--")
        plot(t,yo(t),"b")
        legend("entrada","saida");
        legend('location','northeast');
endfunction

function exemplo_B4o()   %seno Estado Nulo
        ini=-0.2;
        fim=1.0;
        t=linspace(ini,fim,3000);
        xi= @(t) 37*cos(15*t).*f_u(t);
        yo= @(t)(-44.4*exp(-30*t)+18*exp(-90*t)+26.4*cos(15*t)+19.2*sin(15*t)).*f_u(t);
        hold on; grid on;
        plot(t,xi(t),"r");
    %    plot(t,yo(t),"b");
        legend("entrada","saida");
        legend('location','southeast');
endfunction

function exemplo_B4n() %seno Completa
        ini=-0.2;
        fim=1.0;
        t=linspace(ini,fim,3000);
        xi= 0*t;
        yn= @(t)(30*exp(-30*t)-10*exp(-90*t)).*f_u(t);
        hold on; grid on;
        plot(t,xi,"r--")
        plot(t,yn(t),"b")
        legend("entrada","saida");
        legend('location','southeast');
        axis([ini,fim,-40,40]);
endfunction

function exemplo_B4c() %seno Completa
        ini=-0.2;
        fim=1.0;
        t=linspace(ini,fim,3000);
        xi= @(t) 37*cos(15*t).*f_u(t);
        yo= @(t)(-44.4*exp(-30*t)+18*exp(-90*t)+26.4*cos(15*t)+19.2*sin(15*t)).*f_u(t);
        yn= @(t)(30*exp(-30*t)-10*exp(-90*t)).*f_u(t);
        yc= @(t) yn(t)+yo(t);
        hold on; grid on;
        plot(t,xi(t),"r--")
        plot(t,yc(t),"b")
        legend("entrada","saida");
        legend('location','southeast');
endfunction

function y_out=Sistema_Eletrico_2aOrdem(x_in,y1,y2,T)
    den=1+120*T+2700*T^2;
    a1=(2+120*T)/den;
    a2=-1/den;
    b0=(2700*T^2)/den;
    xi=[0 0 x_in];
    yi=[y1 y2];
    for (i=3:length(xi))
         yi(i)=a1*yi(i-1)+a2*yi(i-2)+b0*xi(i);
    end
    y_out=yi(3:length(xi));
endfunction

function exemplo_B5(T)
    ini=0;
    fim=1.0;
    N=ceil((fim-ini)/T +1);
    n=[-1:N-2];
    subplot(521)  %degrau
        x=[0 ones(1,N-1)];
        stem(n*T,x,'r','marker','none');
    subplot(522) %resposta ao degrau
        y=Sistema_Eletrico_2aOrdem(x,0,0,T);
        stem(n*T,y,'b','marker','none');
   subplot(523) %impulso
        x=[0 1 zeros(1,N-2)];
        stem(n*T,x,'r','marker','none');
    subplot(524) %resposta ao impulso
        y=Sistema_Eletrico_2aOrdem(x,0,0,T);
        stem(n*T,y,'b','marker','none');
    subplot(525) %exponencial
        xi = @(t) 10*exp(-15*t).*f_u(t);
        stem(n*T,xi(n*T),'r','marker','none');
    subplot(526) %resposta a exponencial
        y=Sistema_Eletrico_2aOrdem(xi(n*T),0,0,T);
        stem(n*T,y,'b','marker','none');
        axis([-0.2,1,0,10]);
    subplot(527) %cosseno
        xi = @(t) 37*cos(15*t).*f_u(t);
        stem(n*T,xi(n*T),'r','marker','none');
    subplot(528) % resposta ao coseno
        y=Sistema_Eletrico_2aOrdem(xi(n*T),0,0,T);
        stem(n*T,y,'b','marker','none');
    subplot(529) %cosseno
         stem(n*T,xi(n*T),'r','marker','none');
    subplot(5,2,10) % resposta ao coseno
        y=Sistema_Eletrico_2aOrdem(xi(n*T),20,20,T);
        stem(n*T,y,'b','marker','none');
endfunction

function exemplo_C1() %degrau
        ct=1/2;
        ini = -ct;
        fim = 5*ct;
        t=linspace(ini,fim,3000);
        hold on; grid on;
        y_u= @(t) (1-exp(-2*t).*(cos(6*t)-(1/3)*sin(6*t))).*f_u(t);
        plot(t,f_u(t),"r");
        plot(t,y_u(t),"b");
        legend("entrada","saida");
        legend('location','southeast');
        axis([min(t),max(t),-0.2,1.5]);
endfunction

function exemplo_C2()   % impulso
        ct=1/2;
        ini = -ct;
        fim = 5*ct;
        t=linspace(ini,fim,3000);
        h = @(t) exp(-2*t).*(4*cos(6*t)+(16/3)*sin(6*t)).*f_u(t);
        hold on; grid on;
        quiver(0,0,0,4, 'MaxHeadSize', 0.02, 'Color',
           'red', 'LineWidth', 1);
        plot(t,h(t),"b")
        legend("entrada","saida");
        legend('location','southeast');

endfunction

function exemplo_C3o()   % Rampa
        ct=1/2;
        ini = -ct;
        fim = 5*ct;
        t=linspace(ini,fim,3000);
        xi = @(t) t/3.*f_u(t);
        yo = @(t) ((1/3)*t-(1/18)*exp(-2*t).*sin(6*t)).*f_u(t);
        hold on; grid on;
        plot(t,xi(t),"r")
        plot(t,yo(t),"b")
        legend("entrada","saida");
        legend('location','southeast');
endfunction

function exemplo_C4o()   % seno Estado NUlo
        ct=1/2;
        ini = -ct;
        fim = 5*ct;
        t=linspace(ini,fim,3000);
        xi= @(t) 0.377*sin(30*t).*f_u(t);
        yo= @(t) (exp(-2*t).*(0.054*cos(6*t)+0.068*sin(6*t))
                -0.054*cos(30*t)-0.001*sin(30*t)).*f_u(t);
        hold on; grid on;
        plot(t,xi(t),"r")
        plot(t,yo(t),"b")
        legend("entrada","saida");
        legend('location','southeast');
endfunction

function exemplo_C4n()   % seno Entrada Nula
        ct=1/2;
        ini = -ct;
        fim = 5*ct;
        t=linspace(ini,fim,3000);
        yn = @(t) exp(-2*t).*(0.3*cos(6*t)+0.1*sin(6*t)).*f_u(t);
        hold on; grid on;
        plot(t,t*0,"r")
        plot(t,yn(t),"b")
        legend("entrada","saida");
        legend('location','southeast');
        axis([min(t),max(t),-0.4,0.4]);
endfunction

function exemplo_C4c()   % seno Completa
        ct=1/2;
        ini = -ct;
        fim = 5*ct;
        t=linspace(ini,fim,3000);
        xi = @(t) 0.377*sin(30*t).*f_u(t);
        yn = @(t) exp(-2*t).*(0.3*cos(6*t)+0.1*sin(6*t)).*f_u(t);
        yo = @(t) (exp(-2*t).*(0.054*cos(6*t)+0.068*sin(6*t))
             -0.054*cos(30*t)-0.001*sin(30*t)).*f_u(t);
        yc = @(t) yn(t)+yo(t);
        hold on; grid on;
        plot(t,xi(t),"r");
        plot(t,yc(t),"b");
        legend("entrada","saida");;
        legend('location','southeast');;
endfunction

function y_out=Sistema_Mecanico_2aOrdem(x_in,y1,y2,T)
    den=1+4*T+40*T^2;
    a1=(2+4*T)/den;
    a2=-1/den;
    b0=(4*T+40*T^2)/den;
    b1=-4*T/den;
    xi=[0 0 x_in];
    yi=[y1 y2];
    for (i=3:length(xi))
         yi(i)=a1*yi(i-1)+a2*yi(i-2)+b0*xi(i)+b1*xi(i-1);
    end
    y_out=yi(3:length(xi));
endfunction

function exemplo_C5(T)
    ct=1/2;
    ini=0
    fim=5*ct
    N=ceil((fim-ini)/T +1);
    n=[-1:N-2];
    subplot(521)  %degrau
        x=[0 ones(1,N-1)]
        bar(n*T,x,0.5,'r','edgecolor','none');
    subplot(522) %resposta ao impulso
        y=Sistema_Mecanico_2aOrdem(x,0,0,T);
        bar(n*T,y,0.5,'b','edgecolor','none');
    subplot(523) %impulso
        x=[0 1 zeros(1,N-2)]
        bar(n*T,x,0.5,'r','edgecolor','none');
    subplot(524) %resposta ao degrau
        y=Sistema_Mecanico_2aOrdem(x,0,0,T);
        bar(n*T,y,0.5,'b','edgecolor','none');
    subplot(525) %rampa
        xi = @(t) t/3.*f_u(t);
        bar(n*T,xi(n*T),0.5,'r','edgecolor','none');
    subplot(526) %resposta a rampa
        y=Sistema_Mecanico_2aOrdem(x,0,0,T);
        bar(n*T,y,0.5,'b','edgecolor','none');
    subplot(527) %seno
        xi = @(t) 0.377*sin(30*t).*f_u(t);
         bar(n*T,xi(n*T),0.5,'r','edgecolor','none');
    subplot(528) % resposta ao seno
        y=Sistema_Mecanico_2aOrdem(xi(n*T),0,0,T);
        bar(n*T,y,0.5,'b','edgecolor','none');
    subplot(529) %seno
        xi = @(t) 0.377*sin(30*t).*f_u(t);
         bar(n*T,xi(n*T),0.5,'r','edgecolor','none');
    subplot(5,2,10) % resposta ao seno
        y=Sistema_Mecanico_2aOrdem(xi(n*T),0.3,0.3,T);
        bar(n*T,y,0.5,'b','edgecolor','none');
endfunction



