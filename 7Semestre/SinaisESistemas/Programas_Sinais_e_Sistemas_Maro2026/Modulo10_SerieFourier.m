prt=false;

function y=MakePeriodic(f_np,t,t_ini,t_fim)
    tp = mod((t-t_ini),(t_fim-t_ini))+t_ini;
    index = find(tp<t_ini);
    if !isempty(index)
      tp(index)=tp(index) + (t_fim-t_ini);
    end
    y = f_np(tp);
endfunction

function I=intg_simpson(f_trig,t_ini,t_fim,n)
   k=(3000*n)+1; %usa mais pontos para freq. maiores
   x=linspace(t_ini,t_fim,k);
   h=x(2)-x(1);
   y=f_trig(x);
   S=y(1) + 4*sum(y(2:2:k-1)) + 2*sum(y(3:2:k-1)) + y(k);
   I=(h/3)*S;
endfunction

function [F,n]=SerieFourier(f_T,t_ini,t_fim,N)
    T=(t_fim-t_ini);
    w0 = 2* pi  /T;
    C=[];
    for n=1:N
       fsin= @(t) f_T(t).*sin((n-1)*w0*t);
       fcos= @(t) f_T(t).*cos((n-1)*w0*t);
       A=2/T*intg_simpson(fcos,t_ini,t_fim,n);
       B=2/T*intg_simpson(fsin,t_ini,t_fim,n);
       C=[C,(A- i *B)/2];
    end
  %  C=clean(C,1e-10);%zera valores pequenos
    F=[flip(conj(C)) C(2:N)];
    n=[-N+1:N-1];
endfunction

function [A,B,n]=SerieFourier_real(f_T,t_ini,t_fim,N)
    T=(t_fim-t_ini);
    w0 = 2* pi  /T;
    A=[];
    B=[]
    for n=1:N
       fsin= @(t) f_T(t).*sin((n-1)*w0*t);
       fcos= @(t) f_T(t).*cos((n-1)*w0*t);
       A=[A,2/T*intg_simpson(fcos,t_ini,t_fim,n)];
       B=[B,2/T*intg_simpson(fsin,t_ini,t_fim,n)];
    end
    n=[0:N-1];
endfunction

function y=SinteseSerieFourier(F,t,t_ini,t_fim)
    T=(t_fim-t_ini);
    w0 = 2* pi  /T;
    N=(length(F)+1)/2;
    y=F*exp( i *[-N+1:N-1]'*w0*t);
endfunction

function PlotSerieFourier(f,t_ini,t_fim,N)
     clf;
     F=SerieFourier(f,t_ini,t_fim,N);
     T=(t_fim-t_ini);
     w0 = 2* pi  /T;
    subplot(3,1,1)
     t3=linspace(t_ini-2*T,t_fim+2*T,2000);
     grid on; hold on;
     plot(t3,MakePeriodic(f,t3,t_ini,t_fim));
     y3=SinteseSerieFourier(F,t3,t_ini,t_fim); %3 periodos
     plot(t3,y3,'bk');
     t1=linspace(t_ini,t_fim,1000);
     plot(t1,f(t1),'r');  %  1 periodo
    subplot(3,1,2)
     n=[-N+1:N-1];
     stem(n,abs(F),'b','marker','none'); % modulo
     axis([-N,N]);
     grid on;
     title("Modulo"); xlabel("indice"); ylabel("Amplitude");
    subplot(3,1,3)
     fase = atan2(imag(F),real(F));
     index = find((abs(F)<max(abs(F)/1000))|(abs(fase)<1e-4));
     fase(index)=0;
     index = find(fase<=-3);
     fase(index)=fase(index)+2* pi;
     stem(n,fase,'k','marker','none'); % fase
     axis([-N,N]);
     grid on;
     title("Fase");xlabel("indice");ylabel("Fase (rad)");
endfunction

function C = F4(n)
     C = -2*i./(n*pi);
     index=find(mod(n,2)==0);
     C(index)=0;
endfunction;

function C = F5(n)
     C = n*0;
     index=find(mod(n,2)==1);
     C(index)=-20./(n(index).^2*pi^2);
     index=find(n==0);
     C(index)=5;
endfunction;


function F=SerieFourierMediaSinal(y1,fs,fc)
    clf;
    Y1=fft(y1);
    mY1=abs(Y1)/max(abs(Y1));
    fY1=atan2d(imag(Y1),real(y1));
    N=length(y1);
    f0=fs/N;
    delta=floor(fc/f0);
    k=1;
    for indx=1:delta:N/2-delta/2
        ini = floor(indx-delta/2);
        if (ini<1) ini = 1; end
        [mag indx2]=max(mY1(ini:floor(indx-1+delta/2)));
        fase=fY1(indx2+ini-1);
        F(1,k)=(mag*exp( i *fase));
        k=k+1;
        if(k==14) return; end
    end
endfunction

function [F,n,y1,t1]=PlotSerieFourierMediaSinal(y1,fs,fc)
     F=SerieFourierMediaSinal(y1,fs,fc);
     N=length(F);
     F=[conj(flip(F(2:N))) F];
     T=1;
     w0 = 2* pi  /T;
    subplot(2,1,1)
     t3=linspace(-2*T,2*T,2000);
     y3=SinteseSerieFourier(F,t3,-T/2,T/2);
     hold on; grid on;
     plot(t3,y3,'bk');
     t1=linspace(-T/2,T/2,2000);
     y1=SinteseSerieFourier(F,t1,-T/2,T/2);
     plot(t1,y1,'r');
    subplot(2,1,2)
     n=[-N+1:N-1];
     stem(n,abs(F),'b','marker','none'); % modulo
     hold on; grid on;
     title("Modulo");xlabel("indice");ylabel("Amplitude");
endfunction


