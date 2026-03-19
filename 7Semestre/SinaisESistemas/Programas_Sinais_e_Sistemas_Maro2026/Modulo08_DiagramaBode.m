prt = false;

function [M,w]=Resposta_Frequencia_Hs_ilaplace(Hs,w_max)
    Hs=minreal(Hs);
    [num_c,den_c]=tfdata(Hs,'v');
    p=roots(den_c);
    index=find(abs(real(p))>0);
    minp=1;
    if (!isempty(index))  min_p = min(abs(p(index))); end
    ct=1/min_p;
    i=1;
    s=tf('s');
    for (w0=1e-3:w_max/100:w_max)
        Xs=s/(s^2+w0^2);
        ft=ilaplace_fga(Xs*Hs,false);
        t=linspace(10*ct,10*ct+2*pi/w0,1000); %[10ct, 10ct+T];
        M(i)=max(ft(t));
        w(i)=w0;
        i=i+1;
    end
    w =  [w(2:end)-w_max, w];
    M =  [flip(M(2:end)), M];
endfunction

function yp = hval_fga(H, xp)
  [num,den] = tfdata(H,'v');
  N=polyval(num,xp);
  D=polyval(den,xp);
  yp=N./D;
end

w_max=0;
function Resposta_Frequencia_Hs(Hs,w_max,cor)
    if (w_max==0)  w_max=10*max(abs(roots(Hs.den))); end
    w=linspace(-w_max,w_max,1024);
    Hw=hval_fga(Hs,i*w);
   subplot(212)
    plot(w,atan2(imag(Hw),real(Hw)),cor);
    title("Fase"); xlabel("w radianos/s"); ylabel("Fase(rad)");
    axis([-w_max,w_max,-pi,pi]);
    xticks (linspace(-w_max,w_max,9));
    hold on; grid on;
   subplot(211)
    plot(w,(abs(Hw)),cor);
    title("Amplitude"); xlabel("w radianos/s"); ylabel("|H(jw)|");
    axis([-w_max,w_max,min(abs(Hw)),max(abs(Hw))]);
    xticks (linspace(-w_max,w_max,9));
    hold on; grid on;
endfunction

function DiagramaBode(Hs,logf1,logf2,cor)
   logf=linspace(logf1,logf2,1024);
   w=2* pi*(10.^logf);
   Hw=hval_fga(Hs,i*w);
  subplot(211)
   plot(logf,20*log10(abs(Hw)),cor)
   title("Amplitude em dB");
   xlabel("log10(f) Hertz"); ylabel("|H(jw)|dB");
   hold on; grid on;
  subplot(212)
   plot(logf,atan2(imag(Hw),real(Hw)),cor);
   title("Fase");
   xlabel("log10(f) Hertz"); ylabel("Fase(rad)");
   hold on; grid on;
endfunction

function Resposta_Frequencia_Hs4(Hs,logf1,logf2)
    clf;
    subplot(221)
     w=linspace(-2* pi*10.^logf2,2* pi*10.^logf2,1024);
     Hw=hval_fga(Hs, i*w);
     plot(w,abs(Hw));
     axis([min(w)/2,max(w)/2,0,max(abs(Hw ))*1.1]);
     title("Amplitude"); xlabel("w (rad/s)"); ylabel("|H(jw)|");
     grid on;
    subplot(223)
     plot(w,atan2(imag(Hw),real(Hw)));
     title("Fase");xlabel("w (rad/s)"); ylabel("Fase (rad)");
     axis([min(w)/2,max(w)/2,-pi,pi]);
     grid on;
    subplot(222)
     logf=linspace(logf1,logf2,1024);
     w=2* pi*10.^logf;
     Hw=hval_fga(Hs, i*w);
     plot(logf,20*log10(abs(Hw)),'r');
     axis([logf1,logf2,-100,20]);
     title("Amplitude dB"); xlabel("log(f) Hz"); ylabel("|H(jw)|dB");
     grid on;
    subplot(224)
     plot_polos_zeros(Hs);
endfunction

