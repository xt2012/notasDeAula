prt = false;

function y=windows_dsp(N,tipo)
    n=0:N-1;
    if tipo== 1  %retangular
        y=ones(1,N);
    elseif tipo== 2   %bartlett
        y=1-abs((n-(N-1)/2)/((N-1)/2));
    elseif tipo== 3  %hanning
        y=0.5-0.5*cos(2* pi *n/(N-1));
    elseif tipo== 4  %hamming
        y=0.54-0.46*cos(2* pi *n/(N-1));
    elseif tipo== 5  %blackman
        y=0.42-0.50*cos(2* pi *n/(N-1))+0.08*cos(4* pi *n/(N-1));
    elseif tipo== 6  %flattop
        y=(1.0-1.93*cos(2* pi *n/(N-1))+1.29*cos(4* pi *n/(N-1))...
           -0.388*cos(6* pi *n/(N-1))+0.028*cos(8* pi *n/(N-1)))/4;
    elseif tipo== 7  %gausssiana
        y=exp(-0.5*((n-(N-1)/2)/(0.35*(N-1)/2)).^2);
    else tipo== 3  %hanning
        y=0.5-0.5*cos(2* pi *n/(N-1));
    end

endfunction

function windows_dsp_plot()
   N=1024;
   janelas=["retangular";"bartlett";"hanning";"hamming";"blackman";"flat-top";"Gauss"];
   for(k=1:7)
      subplot(7,2,k*2-1);
      y   =windows_dsp(N,k);
      plot(linspace(-0.5,0.5,N+2),[0 y 0],'r');
      axis([-0.5,0.5,-0.1,1.2]);
      title(janelas(k,:));
      grid on;
      subplot(7,2,k*2);
      y   =windows_dsp(N,k);
      w=linspace(1e-12, pi /50,2048);
      Ys=tftd_fga(y,w);
      plot(w/(2* pi ),20*log10(abs(Ys))-20*log10(max(abs(Ys))),'b');
      grid on;
   end
endfunction

function plot_dft_tftd(f,fs,T)
  clf;
   subplot(311) % amostar o sinal f, plotar o y e yn
    yn=sample_signal3(f,fs,0,T,true);
    N=length(yn);
   subplot(312) %calcular e plotar dft de yn
    f_hertz = linspace(-fs/2+mod(N,2)*fs/(2*N), fs/2-fs/N+mod(N,2)*fs/(2*N), N);
    Yn = fft(yn);
    stem(f_hertz,fftshift(abs(Yn)),'k','marker','none');
    title("|Y(m)| = DFT[y(n)]");xlabel("f Hertz");
    grid on;
   subplot(313) %calcular tftd
    w=linspace(-pi,pi,2000);
    hold on; grid on;
    Fs=tftd_fga(yn,w);;
    plot(w/(2*pi)*fs,abs(Fs),'b'); %plot tdft
    stem(f_hertz,fftshift(abs(Yn)),'k','marker','none'); %plot dft
    title("|F(m)| e |Fs(jw)| - 1 periodo da DTFT");xlabel("f Hertz");
endfunction

function ys=Efeito_Janelamento(f,fs,T,janela)
	clf;
    subplot(421) % Plotar sinal continuo de -T a 2T
     ys=sample_signal3(f,fs,0,T,false);
     N=length(ys);
     Nc=3000;
     t=linspace(0,T,Nc);
     yc = f(t);
     grid on; hold on;
     plot([t-T t t+T],[ f([t-T t t+T]) ;[yc* nan  yc yc* nan ]]);
     title("Sinal Continuo yc(t) com extensao [-inf,+inf]");
    subplot(422); % Aproximar a TF do sinal continuo de [-inf,+inf]
     y_inf = f(linspace(-100*T,100*T,200*N));
     N_inf=length(y_inf);
     f_hertz=linspace(-fs/2,fs/2,N_inf);
     plot(f_hertz,fftshift(abs(fft(y_inf))));
     grid on;
     title("Transformada de Fourier de yc(t)");
    subplot(423); %Janela escolhida entre 0 e T
     js=windows_dsp(N,janela);
     grid on; hold on;
     stem(linspace(-T,2*T-T/(3*N),3*N), [js*0 js js*0] ,'marker','none');
     jc=windows_dsp(Nc,janela);
     plot(linspace(-T,2*T-T/(3*Nc),3*Nc), [jc*0 jc jc*0] );
     title("js(t) - Janela no intervalo [0,T] amostrada com fs");;
     axis([-T,2*T,-0.3,max(jc)+0.2]);
    subplot(424); %Plotar TFTD da janela amostrada com fs
     w=linspace(- pi , pi ,100*N); %(1 periodo)
     plot(w/(2* pi )*fs,abs(tftd_fga(js,w)));
     title("TFTD da janela js(t)");
     grid on;
    subplot(425);   % Sinal multiplicado pela janela
     grid on; hold on;
     stem(linspace(-T,2*T-T/(3*N),3*N), [ys*0 ys.*js ys*0] ,'b','marker','none');
     plot(linspace(-T,2*T-T/(3*Nc),3*Nc), [jc*0 yc.*jc jc*0] );
     title("ys(t) - Sinal yc(t) amostrado com fs e multiplicado por js(t)");
    subplot(426); % TFTD do Sinal Janelado
     plot(w/(2* pi )*fs,abs(tftd_fga(ys.*js,w)));
     title("TDFT do sinal ys(t)");
     grid on;
    subplot(427);  % Sinal amostrado com fs na janela 0 a T
     stem(linspace(-T,2*T-T/(3*N),3*N), [ys*0 ys.*js ys*0] ,'b','marker','none');
     title("ys(n)- Sinal yc(t) amostrado com fs e multiplicado por js(t)");
     grid on;;
    subplot(428); % DFT das amostras ys(n)
     f_hertz=linspace(-fs/2,fs/2-fs/N,N);
     stem(f_hertz,fftshift(abs(fft(ys.*js))),'marker','none');
     title("DFT de ys(n)");
     grid on;
endfunction

function Efeito_Janelamento_bat(t1,t2,janela) % para bat.mat
     clf;
     load("./Dados/bat.mat");
     fs=230.4;
     x=bat';
     N=length(x);
     Ts=1/fs;
     T=N*Ts;
     Td=t2-t1;
     no = floor(t1/Ts)+1;
     Nd = floor(Td/Ts);
    clf();
    subplot(321) % Sinal no Dominio do Tempo
     t=[0:Ts:T-Ts];
     plot(t,x);
     title("x(t) - Sinal no Dominio do Tempo");
     grid on;
    subplot(323); % Janela no Dominio do Tempo;
     w=zeros(1,N);
     w1=windows_dsp(Nd,janela);
     w(no:no+Nd-1)=w1(1:Nd);
     plot(t,w);
     title("w(t) - Janela no Dominio do Tempo");
     grid on;
    axis([0,T,min(w)-0.2,max(w)+0.2]);
    subplot(325); % Sinal Janelado no Dominio do Tempo
     xw= x.*w;
     plot(t,xw);
     axis([0,T,min(x),max(x)]);
     title("xw(t) x(t) - Sinal Janelado no Dominio do Tempo");
     grid on;
    subplot(322); %Sinal no Dominio da Frequencia
     X=fft(x);
     f_hertz=[-N/2:N/2-1]'/N*fs;
     plot(f_hertz,fftshift(abs(X)));
     title("X(jw)- Sinal no Dominio da Frequencia");
     grid on;
    subplot(324); %Janela no Dominio da Frequencia
     eixo_w=linspace(- pi /Nd*20, pi /Nd*20,1000);
     W=tftd_fga(w1,eixo_w);
     plot(eixo_w/(2* pi )*fs,abs(W));
     title("W(jw) - Janela no Dominio da Frequencia");
     grid on;
    subplot(326); %Sinal Janelado no Dominio da Frequencia
     f_hertz=[-N/2:N/2-1]'/N*fs;
     XW=fft(xw);
     f_hertz=linspace(-fs/2,fs/2-fs/N,N);
     plot(f_hertz,fftshift(abs(XW)));
     axis([-fs/2,fs/2,0,max(abs(X)/4)]);
     title("X(jw)*W(jw) - Sinal Janelado no Dominio da Frequencia");
     grid on;
endfunction

function Efeito_Janelamento_Mares(t1,t2,janela)
     A=fscanfMat("./Dados/MareNovaYork_1hora.txt");
     fs = 24 % (24 amostras/dia)
     x=A(:,2)';
     x=x-mean(x)
     N=length(x)
     Ts=1/fs
     T=N*Ts;
     Td=t2-t1
     no = int(t1/Ts)+1;
     Nd = int(Td/Ts);
     clf;
    subplot(321)
     t=[0:Ts:T-Ts];
     plot(t,x);
     title("x(t) - Sinal no Dominio do Tempo");
    axis([-0.2,min(x)-0.2;T+0.2,max(x)+0.2]);
    subplot(323);
     w=zeros(1:N)
     w1=windows_dsp(Nd,janela);
     w(no:no+Nd-1)=w1(1:Nd);
     plot(t,w);
    axis([0,min(w)-0.2;T,max(w)+0.2]);
     title("w(t) - Janela no Dominio do Tempo");
     plot([-0.2,t1,t1],[0,0,w(1)])
     plot([t2,t2,T+0.2],[w(Nd),0,0])
    subplot(325);
     xw = x.*w;
     plot(t,xw);
    axis([0,min(x)-0.2;T,max(x)+0.2]);
     title("x(t) w(t) - Sinal janelado no Domino do Tempo");
     plot([-0.2,0,0],[0,0,xw(1)])
     plot([T,T,T+0.2],[x(N),0,0])
    subplot(322);
     X=fft(x);
     f_hertz=[-N/2:N/2-1]'/N*fs
     plot(f_hertz,fftshift(abs(X)));
     title("X(jw) - Sinal no Dominio da Frequencia");
    %axis([-fs/2,0;fs/2,max(abs(X))]);
    subplot(324);
     eixo_w=linspace(- pi /Nd*20, pi /Nd*20,1000)
     W=tftd_fga(w1,eixo_w)
     plot(eixo_w/(2* pi )*fs,abs(W))
    subplot(326);
     f_hertz=[-N/2:N/2-1]'/N*fs
     XW=fft(xw);
     plot(f_hertz,fftshift(abs(XW)));
     title("X(jw)*W(jw) - Sinal Janelado no Dominio da Frequencia");
    axis([-fs/4,0;fs/4,max(abs(X)/16)]);
endfunction

function TF_janelada_octave(y0,fs,Tw,ganho)
     clf;
     Ts=1/fs;
     window = ceil(Tw/Ts)
     pad = ceil(window/2);
     step = ceil(window/20)
     nfft = 2^nextpow2(window);
     [A,f,t]=specgram([zeros(pad,1);y0;zeros(pad,1)], nfft, fs, window, window-step);
     imagesc(t,f,abs(A).^(1/ganho));
     colormap(jet);
     axis xy;
endfunction

function TF_janelada_octave_new(y0, fs, Tw, ganho,tipo)
    clf;
    Ts = 1/fs;
    window = ceil(Tw/Ts);
    pad = ceil(window/2);
    step = ceil(window/20);
    nfft = 2^nextpow2(window);
    w=windows_dsp(window,tipo)';
    [A, f, t]=specgram([zeros(pad,1);y0;zeros(pad,1)], nfft, fs, w, window-step);
    imagesc(t, f, abs(A).^(1/ganho));
    axis xy;
    xlabel('Tempo (s)');
    ylabel('Frequência (Hz)');
    janelas=["retangular";"bartlett";"hanning";"hamming";"blackman";"flat-top";"Gauss"];
    disp(janelas(tipo,:))
    str0 = sprintf('Espectrograma usando janela %s',janelas(tipo,:));
    title(str0);
    colormap(jet);
    colorbar;
endfunction

function TF_janelada(y,fs,Tw,janela,ganho)
     clf;
     N=length(y);
     if (mod(N,2)==1)
       N=N+1;
       y=[y; y(end)];
     endif
     Ts=1/fs;
     T=N*Ts;
     t=0:Tw/20:T;
     f=[0:N/2-1]'/N*fs;
     M=length(t);
     Nd = floor(Tw/Ts);
     tdft = zeros(M,N/2);
     col = 0;
     for(t_ini=0:Tw/20:T)
         col = col+1;
         no = floor((t_ini-Tw/2)/Ts)+1;
         if(no<1) no=1;  end
         n1 = floor((t_ini+Tw/2)/Ts)+1;
         if (n1>N) n1=N; end
         w=zeros(1,N);
         w(no:n1)=windows_dsp(n1-no+1,janela);
         XW=abs(fft(y.*w));
         tdft(col,1:N/2)=(XW(1:N/2));
     end
     colormap(jet);
     imagesc(t,f,abs(tdft').^(1/ganho));
     axis xy;
     shading('interp');
     grid on;
endfunction

function [eixo_f,Ye]=EstimacaoEspectral(y_in,fs)
     y=y_in'-mean(y_in);
     N=length(y);
     passo = floor(N/100);
     if passo<1   passo=1 end
     janela = floor(N/5);
     if janela<4  janela=4 end
     ni=0;
     Ye = zeros(1,N);
     for(i=floor(janela/2)+1:passo:floor(N-janela/2))
         n0 = i;
         n1 = floor(n0+janela/2)+1;
         w=[zeros(1,n0-1) windows_dsp(n1-n0+1,3) zeros(1,N-n1)];
         XW=abs(fft(y.*w));
         Ye=Ye+XW.^2;
         ni=ni+1;
     end
     eixo_f=[-N/2:N/2-1]'*fs/N;
     Ye=sqrt(Ye/ni);
     subplot(311)
         Ts=1/fs;
         t=[0:N-1]*Ts;
         plot(t,y);
         grid on;
         axis([min(t),max(t)]);
         title("Sinal no Tempo");
     subplot(312)
         plot(eixo_f,fftshift(Ye)/max(Ye));
         grid on;
         title("Estimacao Espectral");
     subplot(313)
         Ys=abs(fft(y));;
         plot(eixo_f,fftshift(Ys)/max(Ys));
         grid on;
         title("DFT");
endfunction

function [eixo_f,R]=EstimacaoCovariancia(x,fs,janela)
   N=length(x)
   x=x-mean(x)
   for (m=0:N-1)
     x2=circshift(x,m)
     r(m+1)=x*x2'
   end
   r=r/N;
   w=windows_dsp(N,janela)';
   R=sqrt(fftshift(abs(fft(r.*w))))
   eixo_f=[-N/2:N/2-1]'*fs/N;
   plot(eixo_f,R)
endfunction

