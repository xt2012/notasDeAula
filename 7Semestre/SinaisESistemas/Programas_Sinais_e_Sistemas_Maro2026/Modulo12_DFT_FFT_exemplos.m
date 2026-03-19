prt = false;

function anima1()
    figure(1);
    clf
    dt=0.025
    i=1
    for t=0:dt:8.8-0.5
        Efeito_Janelamento_bat(t,t+0.5,3)
        nome=sprintf("test%d",i)
        xs2jpg(0,nome)
        i=i+1
    end
    animaGIF("test*.jpeg","fft_bat.gif",50)
endfunction

function anima2()
    figure(1);
    clf;
    dt=5
    i=1
    for t=0:dt:350-50
        Efeito_Janelamento_Mares(t,t+50,3)
        nome=sprintf("test%d",i)
        xs2jpg(0,nome)
        i=i+1
    end
    animaGIF("test*.jpeg","fft_mare.gif",70)
endfunction

function exemplo_dft_temperatura()
    [t,y]=textread("./Dados/temperatura_DF_2.txt","%f %f",686);
    y0=y'-mean(y);
    N=length(y);
    fs = 12; % (12 amostras/ano)
  subplot(411)
    plot(t,y0);
    axis([min(t),max(t)]);
    grid on;
    title("Temperatura mensal do DF");
  subplot(412)
    w=linspace(-pi , pi ,3600);
    Ys=tftd_fga(y0,w);
    plot(w/(2* pi )*fs,abs(Ys));
    axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
    grid on;
    title("TDFT");
  subplot(413)
    f=linspace(-fs/2,fs/2,N);
    Y0=dft_fga(y0);
    stem(f,fftshift(abs(Y0)),'b','marker','none');
    axis([min(f),max(f)]);
    grid on;
    title("DFT");
  subplot(414)
    y_inverso=idft_fga(Y0);
    plot(t,y_inverso+mean(y),'r');
    grid on;
    axis([min(t),max(t)]);
    title("IDFT");
endfunction

function exemplo_dft_morcego()
    load("./Dados/bat.mat");
    y=bat';
    N=length(y);
    fs = 230.4; % (230.4khz)
    Ts = 1/fs;
    t=[0:Ts:N*Ts-Ts];
  subplot(411)
    plot(t,y);
    grid on;
	axis([min(t),max(t)]);
	title("Som Morcego");
  subplot(412)
    w=linspace(-pi,pi ,3600);
    Ys=tftd_fga(y,w);
    plot(w/(2* pi )*fs,abs(Ys));
    grid on;
   	axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
	  title("TDFT");
  subplot(413)
    f=linspace(-fs/2,fs/2,N);
    Y=dft_fga(y);
    stem(f,fftshift(abs(Y)),'b','marker','none');
    grid on;
  	axis([min(f),max(f)]);
	  title("DFT");
  subplot(414)
    y_inverso=idft_fga(Y);
    plot(t,y_inverso,'r');
    grid on
	  axis([min(t),max(t)]);
	  title("IDFT");
endfunction

function exemplo_dft_MareNovaYork_1h()
    A=load("./Dados/MareNovaYork_1hora.txt");
    y=A(:,2);
    n=A(:,1);
    N=length(y);
    fs = 24; % (24 amostras/dia)
    f=linspace(-fs/2,fs/2-fs/N,N);
    subplot(311)
    plot(n,y);
    title("Mare Nova York");
    xlabel("dia");
    ylabel("metros");
    fprintf("cheguei aqui")
    hamm=windows_dsp(N,4)';
    Y=fft((y-mean(y)).*hamm);
    subplot(312)
    stem(f,fftshift(abs(Y)),'b','marker','none');
    title("Analise Frequencia Mare Nova York");
    xlabel("Freq=Oscilacoes//dia");
    ylabel("|T(jw)|");
    w=linspace(-pi , pi ,3600);
	  axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
    subplot(313)
    stem(f,fftshift(abs(Y)),'b','marker','none');
    axis([-fs/256,fs/256,-20,max(abs(Y))/8])
    title("Analise Frequencia Mare Nova York");
    xlabel("Freq=Oscilacoes/dia");
    ylabel("|T(jw)|");
endfunction

function exemplo_dft_MareNovaYork()
    A=load("./Dados/MareNovaYork.txt");
    y=A(:,2);
    n=A(:,1);
    N=length(y);
    fs = 12; % (12 amostras/ano)
    f=linspace(-fs/2,fs/2-fs/N,N);
    subplot(311)
    plot(n,y)
    title("Mare Nova York");
    xlabel("ano");
    ylabel("metros");
    hamm=windows_dsp(N,4)';
    Y=fft((y-mean(y)).*hamm);
    subplot(312)
    stem(f,fftshift(abs(Y)),'b','marker','none');
    w=linspace(-pi , pi ,3600);
	  axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
    title("Analise Frequencia Nova York");
    xlabel("Freq=Oscilacoes/ano");
    ylabel("|T(jw)|");
    subplot(313)
    stem(f+fs/(2*N),fftshift(abs(Y)),'b','marker','none');
    axis([-0.1,0.1,-20,max(abs(Y))/4])
    title("Analise Frequencia Nova York");
    xlabel("Freq=Oscilacoes/ano");
    ylabel("|T(jw)|");
endfunction

function exemplo_dft_EuroDolar()
    A=load("./Dados/Euro-dolar.txt");
    y=A(:,2);
    n=A(:,1);
    N=length(y);
    fs = 365; % (365 amostras/ano)
    f=linspace(-fs/2,fs/2-fs/N,N);
    subplot(311)
    plot(n,y)
    title("Cotacao do Euro em Dolar");
    xlabel("ano");
    ylabel("Dolar");
    hamm=windows_dsp(N,4)';
    Y=fft((y-mean(y)).*hamm);
    subplot(312)
    stem(f,fftshift(abs(Y)),'b','marker','none');
    title("Analise Frequencia Cotacao do Dolar");
    xlabel("Freq=Oscilacoes/ano");
    ylabel("|Y(jw)|");
    w=linspace(-pi , pi ,3600);
	  axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
    subplot(313)
    stem(f,fftshift(abs(Y)),'b','marker','none');
    axis([-2.5,2,5,-20,max(abs(Y))])
    title("Analise Frequencia Cotacao do Dolar");
    xlabel("Freq=Oscilacoes/ano");
    ylabel("|Y(jw)|");
endfunction

function exemplo_dft_golfinho()
    [y,fs,bits]=wavread("./dados/dol17.wav");
    N=length(y)
    fs = fs/1000  % (kHz)
    Ts = 1/fs;
    n=[0:Ts:N*Ts-Ts];
    f=linspace(-fs/2,fs/2-fs/N,N);
    subplot(211)
    plot(n,y)
    title("Som Golfinho");
    xlabel("ms");
    ylabel("Amplitude");
    Y=fft(y-mean(y));
    subplot(212)
    stem(f,fftshift(abs(Y)),'b','marker','none');
	axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
    title("Analise Frequencia Golfinho");
    xlabel("Freq em kHz");
    ylabel("|T(jw)|");
endfunction


function exemplo_dft_baleia()
    [y,fs,bits]=wavread("./dados/WHALEWHI.wav");
    N=length(y) ;
    fs = fs/1000 ;  % (kHz)
    Ts = 1/fs;
    n=[0:Ts:N*Ts-Ts];
    f=linspace(-fs/2,fs/2-fs/N,N);
    subplot(211);
    plot(n,y);
    title("Assovio Baleia");
    xlabel("ms");
    ylabel("Amplitude");
    hamm=windows_dsp2(N,4);
    Y=fft((y-mean(y)).*hamm);
    subplot(212);
    stem(f,fftshift(abs(Y)),'b','marker','none');
	axis([min(w/(2* pi )*fs),max(w/(2* pi )*fs)]);
    title("Analise Frequencia Assovio Baleia","Freq em kHz","|T(jw)|");
endfunction

function exemplo_dft_music()
    [y0,fs]=audioread("./dados/music2.wav");
    N=length(y0);
    Ts = 1/fs;
    t=[0:Ts:N*Ts-Ts];
  subplot(311)
    plot(t,y0);
    grid on;
	  axis([min(t),max(t)]);
	  title("Musica");
  subplot(312)
    f=linspace(-fs/2,fs/2,N);
    Y=fft(y0);
    stem(f,fftshift(abs(Y)),'b','marker','none');
    axis([-fs/2,fs/2-fs/N,0,max(abs(Y))])
	title("Analise Frequencia Musica");
    grid on
  subplot(313)
    plot(f,fftshift(abs(Y)),'r');
    axis([-500,500,0,max(abs(Y))])
	title("Analise Frequencia Musica");
    grid on
endfunction

function speed_test()
    tp1(1:13)=tp2(1:13)=tp3(1:13)=tp4(1:13)=tp5(1:13)=0;
    N=100;
    for (n=1:N)
		printf("\n%d\n,",n);
        for k=4:13
		    printf("%d ",k);
            x=linspace(0,1,2^k);
            tic();  dft_fga(x)    ;dt=toc();tp1(k)=tp1(k)+dt;
            tic();  fft_rec_fga(x);dt=toc();tp2(k)=tp2(k)+dt;
            tic();  fft_dif_fga(x);dt=toc();tp3(k)=tp3(k)+dt;
            tic();  fft_dit_fga(x);dt=toc();tp4(k)=tp4(k)+dt;
            tic();  fft(x)        ;dt=toc();tp5(k)=tp5(k)+dt;
        end
    end
    k=[4:13];
    hold on; grid on;
    plot(k,log2(tp1(k)/N),'b');
    plot(k,log2(tp2(k)/N),'y');
    plot(k,log2(tp3(k)/N),'g');
    plot(k,log2(tp4(k)/N),'r');
    plot(k,log2(tp5(k)/N),'bk');
    legend("dft_fga","fft_rec_fga","fft_dif_fga","fft_dit_fga",
	         "fft octave",'location','northwest');
	title("FFT Speed Test");
	xlabel("10^n amostras");ylabel("tempo 1/10^n segundos");
	axis([4,13]);xticks([1:13]);yticks([-15:5])
endfunction


