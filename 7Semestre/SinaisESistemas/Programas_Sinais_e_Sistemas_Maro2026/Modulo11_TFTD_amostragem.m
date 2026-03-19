prt=false;

function [yn, n]=sample_signal(f,fs,T,prt)
    [yn, n]=sample_signal3(f,fs,0,T,prt);
endfunction

function [yn, n]=sample_signal3(f,fs,t1,t2,prt)
    Ts=1/fs;
    P=t2-t1;
    N = floor(P/Ts);
    if mod(N,2) == 0   N=N+1 end
    n=[0:N-1];
    ts = n*Ts + t1;
    yn = f(ts);
    if (prt)
        printf("N=%d n=[0:%d] fs=%f Ts=%f P=%f\n",N,N-1,fs,Ts,P);
        t=linspace(t1,t2,2000);
        yt = f(t);;
        grid on; hold on;
        stem(ts,yn,'k','marker','none');
        plot(t,yt);
        axis([t1,t2]);
    end
endfunction

function [yn, n]=sample_signal2(f,fs,t1,t2,prt)
    Ts=1/fs;
    P=t2-t1;
    N = floor(P/Ts)
    if mod(N,2) == 1   N=N+1; end
    n=[0:N-1];
    ts = n*Ts + t1;
    yn = f(ts);
    if (prt)
        printf("N=%d n=[0:%d] fs=%f Ts=%f P=%f\n",N,N-1,fs,Ts,P);
        t=linspace(t1,t2,2000);
        yt = f(t);;
        grid on; hold on;
        stem(ts,yn,'b','marker','none');
        plot(t,yt);
        axis([t1,t2]);
    end
endfunction

function Ys=tftd_fga(y,w)
     N=length(y);
     n=[0:N-1];
     Ys=y*exp(- i *n'*w);
endfunction

function Ys=tftd_fga2(y)
     Ys=y*exp(- i *[0:length(y)-1]'*linspace(-pi,pi,3600));
endfunction

function Resposta_Frequencia_yn(yn,ws,cor,nper)
   w=linspace(-nper*pi,nper*pi,nper*1024);
   H =tftd_fga(yn,w);
   %spi = part(%chars.greek.lower, 19);
   spi = "pi";
   plot(w/(2*pi)*ws,abs(H)/max(abs(H)),cor);
   axis=[-ws/2,ws/2];
   ax1=gca;
   set(ax1, 'XAxisLocation', 'top','YAxisLocation','Right');
   ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
   xticks (ax2,[-pi:pi/4:pi]/(2*pi)*ws);
   set(ax2, 'XAxisLocation', 'bottom','YAxisLocation','Left');
   set(ax2, 'XLim', get(ax1, 'XLim'));
   set(ax2, 'XTickLabel', {'-pi','-3pi/4','-pi/2','-pi/4','0','pi/4','pi/2','3pi/4','pi'});
   title("|H(jw)|");
   grid on;
endfunction

function ys=itftd_fga(Y,w)
   N=length(Y);
   M=length(w);
   dw = (max(w)-min(w))/M;
   for n=0:N-1;
       ys(n+1)=0;
       for m=0:M-1;
           ys(n+1)=ys(n+1)+Y(m+1)*exp( i *n*w(m+1))*dw;;
       end
    end
    ys=1/(2* pi  )*clean(ys);;
endfunction


function Ys=tftd_fga_n0(y,n0,w)
     N=length(y);
     n=[n0:N+n0-1];
     Ys=y*exp(- i *n'*w);
endfunction

function Ys=PlotTFTD_amostras2(ys,n0)
   subplot(311)
    n=[n0:length(ys)+n0-1];
  	stem(n,ys,'b','marker','none');
    grid on;
	title("Sinal Amostrado no Tempo");
   subplot(312)
    w=linspace(- pi  , pi  ,5000);
    Ys=tftd_fga_n0(ys,n0,w);
    plot(w,abs(Ys)/ max(abs(Ys)));
    grid on;
	title("Modulo TFTD");
    axis([- pi,pi,-0.1,1.1]);
   subplot(313)
    plot(w,atan2(clean(imag(Ys)),clean(real(Ys))));
	title("Fase TFTD");
    grid on;
end

