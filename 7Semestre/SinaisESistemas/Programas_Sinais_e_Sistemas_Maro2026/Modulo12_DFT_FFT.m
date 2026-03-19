prt = false;

function dft_circulo(N)
   clf;
   W=exp(- i *2* pi  *[0:N-1]/N)';
   hold on; grid on;
   title(strcat('Circulo DFT N=',num2str(N)), "fontsize", 20);
   plot(real(W),imag(W),'o');
   angle = linspace(0, 2* pi , 1024);
   plot(cos(angle),sin(angle),'r--');
   for (j=1:N)
     text(real(W(j))+0.05,imag(W(j)),num2str(cleancomplex(W(j)),4), "fontsize", 20)
     plot([0,real(W(j))],[0,imag(W(j))],'b--');
   endfor
   axis ([-1.2,1.2,-1.2,1.2]);
endfunction

function W=dft_matrix(N)
    m=0:N-1;
    n=0:N-1;
    W=exp(- i *2* pi  *n'*m/N);
  % W=clean(W);
endfunction

function W=plot_dft_matrix(N)
    for i=1:N
        delta = 360/N * (i-1)
        for j=1:N
            subplot(N,N,(j-1)*N+i)
            grid on
            angulo = -delta*(j-1)
            t = linspace(0,360,1000)
            plot(cosd(t),sind(t))
            plot([-1 1],[0 0])
            plot([0 0],[-1 1])
            gca().axes_visible=["off","off","off"]
            xarrows([0,cosd(angulo)],[0,sind(angulo)],10,5)
        end
    end
    format('v',6)
    W=dft_matrix(N)
    disp(clean(W))
    format('v',12)
endfunction

function Y=dft_fga(y)
    N=length(y);
    W=dft_matrix(N);
    Y=y*W;;
endfunction

function y=idft_fga(Y)
    N=length(Y);
    W=dft_matrix(N);
    y=1/N*Y*conj(W);
endfunction

function X = fft_rec4_fga(x)
    N = length(x);
    if N == 4
        X(1,1) = (x(1)+x(3)) +   (x(2)+x(4));
        X(1,2) = (x(1)-x(3)) - i *(x(2)-x(4));
        X(1,3) = (x(1)+x(3)) -   (x(2)+x(4));
        X(1,4) = (x(1)-x(3)) + i *(x(2)-x(4));
    else
        N2 = N/2;
        X_par   = fft_rec_fga(x(1:2:N-1));
        X_impar = fft_rec_fga(x(2:2:N));
        W = exp(-2 *  pi   *  i  *(0:N2-1)/N) ;
        X_impar = W.* X_impar;
        X = [ X_par + X_impar , X_par - X_impar ];
    end
endfunction

function X = fft_rec_fga(x)
    N = length(x);
    if N == 1
        X = x;
    else
        N2 = N/2;
        X_par   = fft_rec_fga(x(1:2:N-1));
        X_impar = fft_rec_fga(x(2:2:N));
        W = exp(-2 *  pi   *  i  *(0:N2-1)/N) ;
        X_impar = W.* X_impar;
        X = [ X_par + X_impar , X_par - X_impar ];
    end
endfunction

function x = ifft_rec_fga(X)
    N = length(X);
    x=1/N*fft_rec_fga(conj(X));
endfunction

function X = fft_dit_fga(x)
    M=ceil(log2(length(x))); %numero estagios
    x=[x zeros(1,(2^M)-length(x))];
    N=length(x);
    N2=1;
    X=bitrevorder(x);
    for estagio=1:M;  %estagios
        n=0:(N2-1);
        W=exp(- i *2* pi  *(2^(M-estagio))*n/N);
        for index=0:(2^estagio):(N-1); %borboletas
             k=n+index+1;
             u=X(k);
             v=X(k+N2).*W;
             X(k)   = u+v;
             X(k+N2)= u-v;
        end;
        N2=2*N2;
    end;
endfunction

function x = ifft_dit_fga(X)
    N = length(X);
    x=1/N*fft_dit_fga(conj(X));
endfunction

function X = fft_dif_fga(x)
    M=ceil(log2(length(x))); %numero estagios
    x=[x zeros(1,(2^M)-length(x))];
    N=length(x);
    N2=N/2;
    for estagio=1:M;  %estagios
        n=0:(N2-1);
        W=exp(- i *2* pi  *(2^(estagio-1))*n/N);
        for index=0:(N/(2^(estagio-1))):(N-1); %borboletas
                k=n+index+1;
                u= x(k)+x(k+N2);
                v=(x(k)-x(k+N2)).*W;
                x(k)   =u;
                x(k+N2)=v;
        end;
        N2=N2/2;
    end;
    X=bitrevorder(x);
endfunction

function x = ifft_dif_fga(X)
    N = length(X);
    x=1/N*fft_dif_fga(conj(X));
endfunction

function PlotDFT_amostras(yn)
    subplot(211)
     stem([0:length(yn)-1],yn,'r','marker','none');
     title("f(n) Discreto e nao-Periodico");;
     xlabel("n"); ylabel("f(n)");
    subplot(212)
     Ys=fftshift(fft(yn));
     DFT=abs(Ys)/ max(abs(Ys));
     eixo_w = linspace(- pi  , pi  ,length(Ys));
     plot(eixo_w,DFT);
     title("Sinal no Dominio da Frequencia");
	   xlabel("w normalizada");
     ylabel(" % |F(jw)| ");
     grid on
end

function PlotDFT_amostras2(yn,fs)
    subplot(211)
     plot([0:length(yn)-1],yn,0.2)
     title("f(n) Discreto e nao-Periodico");
     xlabel("n");
     ylabel("f(n)");
    subplot(212)
     Ys=fftshift(fft(yn))
     DFT=abs(Ys)/ max(abs(Ys))
     eixo_w = linspace(- pi  , pi  ,length(Ys))
     plot(eixo_w/(2* pi  )*fs,DFT)
     title("Sinal no Dominio da Frequencia");
	 xlabel(" % |F(jw)| ");
     grid on
end

function yf=System_lp(y,f0,fs)
    N=length(y);
    Y=fft(y);;
	eixo_f=linspace(-fs/2,fs/2, N);
    H=Y*0;
    index=find(abs(eixo_f)<=f0);
    H(index)=1;
	%plot(eixo_f,H);
    Yf = Y .* fftshift(H);
	%plot(eixo_f,fftshift(abs(Yf)));
	yf=real(ifft(Yf));
endfunction

function yf=System_hp(y,f0,fs)
    N=length(y);
    Y=fft(y);;
	eixo_f=linspace(-fs/2,fs/2, N);
    H=Y*0;;
    index=find(abs(eixo_f)>=f0);
    H(index)=1;
	%plot(eixo_f,H);
    Yf = Y .* fftshift(H);
    %plot(eixo_f,fftshift(abs(Yf)));
    yf=real(ifft(Yf));
endfunction

function yf=System_bp(y,f0,f1,fs)
    N=length(y);
    Y=fft(y);
	eixo_f=linspace(-fs/2,fs/2, N);
    H=Y*0;
    index=find((abs(eixo_f)>=f0)&(abs(eixo_f)<=f1));
    H(index)=1;
    %plot(eixo_f,H);
    Yf = Y .* fftshift(H);
	%plot(eixo_f,fftshift(abs(Yf)));
    yf=real(ifft(Yf));
endfunction

function yf=System_br(y,f0,f1,fs)
    N=length(y);
    Y=fft(y);;
	eixo_f=linspace(-fs/2,fs/2, N);
    H=Y*0;
    index=find((abs(eixo_f)<=f0)|(abs(eixo_f)>=f1));
    H(index)=1;
    %plot(eixo_f,H);
    Yf = Y .* fftshift(H);
    %plot(eixo_f,fftshift(abs(Yf)));
    yf=real(ifft(Yf));
endfunction

function h=sinc_filter_fga(fd,N)
    n=[-N/2:N/2-1];
    h=sinc(2*fd*n);
    % h=sinc(2* pi *fd*n).* windows_dsp(K,4)
    h=h/sum(h);
end

function yf=System_conv_lp(y,f0,fs)
    fd=f0/fs;
    K=floor(2/0.01); %passagem 1% de fs
    h_lp=sinc_filter_fga(fd,2*K);
    yf=conv(y,h_lp,'same');
endfunction

%filtro complementar
function yf=System_conv_hp(y,f0,fs)
    fd=f0/fs;
    K= floor(2/0.01); %passagem 1% de fs
    h_hp=-sinc_filter_fga(fd,2*K);
    h_hp(K+1)=h_hp(K+1)+1;
    yf=conv(y,h_hp,'same');
endfunction

% inversao espectral
function yf=System_conv_hp2(y,f0,fs);
    fd=(fs/2-f0)/fs; %passagem 1% de fs
    K=floor(2/0.01);
    n=[-K:K-1];
    h_hp= sinc_filter_fga(fd,2*K).*(-1).^(n);
    yf=conv(y,h_hp,'same');
endfunction

function yf=System_conv_bp(y,f0,f1,fs)
    fd1=f1/fs;
    K=floor(2/0.005); % passagem 0.5% de fs
    h_lp1=sinc_filter_fga(fd1,2*K);
    fd0=f0/fs;
    h_lp0=sinc_filter_fga(fd0,2*K);
    h_bp = h_lp1 - h_lp0;
    yf=conv(y,h_bp,'same');
endfunction

function yf=System_conv_br(y,f0,f1,fs)
    fd1=(fs/2-f1)/fs;
    K=floor(2/0.005);% passagem 1% de fs
    n=[-K:K-1];
    h_hp1= sinc_filter_fga(fd1,2*K).*(-1).^(n);
    fd0=f0/fs;
    h_lp0=sinc_filter_fga(fd0,2*K);
    h_br = h_hp1 + h_lp0;
    yf=conv(y,h_br,'same');
endfunction

f1=0
function Plot_system(tipo,y,fs,f0,f1)
    %tipo 'lp','hp','bp','br'
    %f1 so e usada em 'bp' e 'rp
    Ts = 1/fs;
    N=length(y);
  subplot(321)
    str1=strcat("Sinal na Frequencia " ,num2str(N)," amostras");
	f=linspace(-fs/2+mod(N,2)*fs/(2*N),
         	fs/2-fs/N+mod(N,2)*fs/(2*N), N);
    Y=fft(y);
    stem(f,abs(fftshift(Y)),'b','marker','none');
	title(str1);
	grid on;
  subplot(322)
    str2=strcat("Sinal no Tempo ",num2str(N)," amostras");
    n=[0:N-1]*Ts;
    plot(n,y);
	  title(str2);
	  grid on;
  subplot(323)
    title("Filtro LP/HP/BP/RB na Frequencia")
    if     (tipo=='lp')
         index=find(abs(f)<=f0);  zoom=4/0.01;
    elseif (tipo=='hp')
         index=find(abs(f)>=f0);  zoom=4/0.01;
    elseif (tipo=='bp')
         index=find((abs(f)>=f0)&(abs(f)<=f1)); zoom=4/0.005;
    elseif (tipo=='br')
         index=find((abs(f)<=f0)|(abs(f)>=f1)); zoom=4/0.02;
    else index = []
    end
    H=Y*0;% notar que H esta com shift
    H(index)=1;
    plot(f,abs(H));
    axis([min(f),max(f),-0.5,2]);
    grid on;
    title("Sinal Filtrado No Tempo")
  subplot(324)
    str3=strcat("h[n] - zoom com ",num2str(floor(zoom))," amostras");
    h=real(ifft(fftshift(H)));
    hold on; grid on;
    stem(n,fftshift(h),'b','marker','none');
    plot(n,fftshift(h),'b');
	  title(str3);
    axis([N*Ts/2-zoom*Ts/2,N*Ts/2+zoom*Ts/2,min(h),max(h)]);
  subplot(325)
    Yf = Y .* fftshift(H);
    stem(f,abs(fftshift(Yf)),'b','marker','none');
	  title("Sinal Filtrado Na Frequencia");
	  grid on;
  subplot(326)
    yf=real(ifft(Yf));
    plot(n,yf);
	  title("Sinal Filtrado No Tempo");
	  grid on;
endfunction


