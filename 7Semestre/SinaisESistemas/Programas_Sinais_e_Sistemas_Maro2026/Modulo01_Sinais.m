prt=false;

% Degrau unitario
function  x=f_u(t)
   x =double(t>=0);
end;

% Signum
function  x=f_s(t)
   x =2*f_u(t)-1;
end;

% Funcao Rampa
function  x=f_r(t)
   x =t.*f_u(t);
end;

%pulso
function  x=f_p(t,T)
   x =f_u(t+T/2)-f_u(t-T/2);
end;

%Pulso Area Unitaria largura T
function  x=f_pd(t,T)
   x =1/T*f_p(t,T);
end;

% Impulso Discreto
function  x=f_d(n)
   x = double(n==0);
end;

%Pulso triangular largura T
function  x=f_pt(t,T)
   x =(1-2/T*abs(t)).*f_p(t,T);
end;

% sinc sinais e sistemas
function  x=f_sc(t)
   x =sinc(t);
end;

% Pulso Sinc Periodo T
function  x=f_psc(t,T)
   x =1/T*sinc(t/T);
end;

 % Pulso Gaussiano (med,std)
function  x=f_pg(t,med_in,std_in)
   x =1/(std_in*sqrt(2*pi))*exp(-1/2*((t-med_in)/std_in).^2 );
end;

 % Sinal Exponencial
function  x=f_exp(t,a)
   x = exp(a*t);
end;

 % Pulso Exponencial
function  x=f_pexp(t,T)
   x =1/(2*T)*exp(-abs(t)/T);
end;

function plot_arrow(xo,yo,c,cor)
  u = c/10 * cos(pi/6);
  v = c/10 * sin(pi/6);
  line([xo,xo,xo-u,xo,xo+u], [yo,yo+c,yo+c-v,yo+c,yo+c-v], 'Color', cor, 'LineWidth', 2);
  end

function y = Si(x)
    i_l = find(abs(x)<21.5)
    N = 35 % aprox x pequeno
    c_pol = (-1).^(0:N)./(1:2:2*N+1)^2./factorial(0:2:2*N);
    p_pol=horner(poly(c_pol, "s", "c"), x(i_l).^2)
    y(i_l) = x(i_l).*p_pol;

    i_h = find(abs(x)>=21.5)
    M = 21 % aprox x grande
    sig = (-1).^(0:M);
    c_cos = factorial(0:2:2*M).*sig;
    c_sen = factorial(1:2:2*M+1).*sig;
    p_cos = horner(poly(c_cos,"s","c"), 1./x(i_h).^2)
    p_sen= horner(poly(c_sen,"s","c"), 1./x(i_h).^2)
    y(i_h) = sign(x(i_h))*pi/2 -p_cos.*cos(x(i_h))./x(i_h)...
             - p_sen.*sin(x(i_h))./x(i_h).^2
endfunction

function  y=chirp(fs)
    T=1/fs;
    t=[0:T:10];
    y= sin(pi*t^2) .* (0.5-0.5*cos(2*pi*t/10));
endfunction

function plot_sinal_complexo2(t,y,cont)
     subplot(221)
	     hold on; grid on;
       if (cont)   plot(t,real(y))
       else   plot2d3(t,real(y)) end
       title('real(x)')
     subplot(223)
	     hold on; grid on;
       if (cont)  plot(t,imag(y));
       else  plot2d3(t,imag(y)) end
       title('imag(x)')
     subplot(222)
	     hold on; grid on;
       if (cont)  plot(t,abs(y));
       else  plot2d3(t,abs(y)) end
       title('abs(x)')
       lim=max(abs(y))
       axis( [min(t),max(t),-lim,lim] )
     subplot(224)
	      hold on; grid on;
       if (cont)   plot(t,atan(imag(y),real(y)));
       else     plot2d3(t,atan(imag(y),real(y)));  end
       title('fase(x)')
       axis = [min(t),max(t),-2*pi,2*pi]
endfunction

function plot_sinal_complexo(f,t)
     subplot(221)
       hold on; grid on;
       plot(t,real(f(t)));
       title('real(x)')
     subplot(223)
       hold on; grid on;
       plot(t,imag(f(t)));
       title('imag(x)')
     subplot(222)
       hold on; grid on;
       plot(t,abs(f(t)));
       title('abs(x)')
       lim=max(abs(f(t)));
       axis( [min(t),max(t),-0.1*lim,1.1*lim] )
     subplot(224)
       hold on; grid on;
       plot(t,atan2(imag(f(t)),real(f(t))));
       title('fase(x)')
       axis ([min(t),max(t),-2*pi,2*pi])
endfunction

function plot_sinal_complexo_d(f,t)
     subplot(221)
       hold on; grid on;
       bar(t,real(f(t)),0.3,'b','edgecolor','none');
       title('real(x)')
       axis([min(t),max(t)]);
     subplot(223)
       hold on; grid on;
       bar(t,imag(f(t)),0.3,'b','edgecolor','none');
       title('imag(x)');
       axis([min(t),max(t)]);
     subplot(222)
       hold on; grid on;
       bar(t,abs(f(t)),0.3,'b','edgecolor','none');
       title('abs(x)')
       lim=max(abs(f(t)))*1.2;
       axis( [min(t),max(t),-lim,lim] )
     subplot(224)
       hold on; grid on;
       bar(t,atan2(imag(f(t)),real(f(t))),0.3,'b','edgecolor','none');
       title('fase(x)')
       axis ([min(t),max(t),-2*pi,2*pi])
endfunction

 %fasores(120,45)
 function fasores(N)
    z = [1:N]
    theta = ((z-1)/N)*720;
    xi = zeros(1,N);
    xf = cosd(theta);
    yi = zeros(1,N);
    yf = sind(theta);
    plot3d([-1 1 1 1],[-1 -1 -1 1],[-1 -1 2 2])
    xarrows([xi;xf],[yi;yf],[z;z],0.7, 32*rand(1,N))
 endfunction

 function  plot_arrow_fga(xi ,yi,cor)
  ax=axis;
 % w=get(axes(),'Position')
  w= [0.1300   0.1100   0.7750   0.8150];
  xo = [(xi(1)-ax(1))/(ax(2)-ax(1))*w(3)+w(1),(xi(2)-ax(1))/(ax(2)-ax(1))*w(3)+w(1)];
  yo = [(yi(1)-ax(3))/(ax(4)-ax(3))*w(4)+w(2),(yi(2)-ax(3))/(ax(4)-ax(3))*w(4)+w(2)];
  annotation ("arrow",xo,yo,'color',cor);
 endfunction

  function  plot_fasor(z,h)
    x = real(z);
    y = imag(z);
    colors = ['r', 'g', 'b', 'm','c','y','bk'];
    N= length(z);
    figure;
    hold on;
    grid on;
    arrowScale = 0.1;
    for k = 1:N
      plot([0, x(k)], [0, y(k)], colors(mod(k,6)), 'LineWidth', 1.5);
    end

    for k = 1:N
      theta = atan2(y(k), x(k));  % Direction of the arrow
      ahX = [x(k), x(k) - h*cos(theta + pi/6), x(k) - h*cos(theta - pi/6)];
      ahY = [y(k), y(k) - h*sin(theta + pi/6), y(k) - h*sin(theta - pi/6)];
      fill(ahX, ahY, colors(k));
    end

   % axis equal;
   axis([min([x,y]),max([x,y]),min([x,y]),max([x,y])])
 endfunction

  function  plot_darrow_fga(xi ,yi,cor)
  ax=axis;
 % w=get(axes(),'Position')
  w= [0.1300   0.1100   0.7750   0.8150];
  xo = [(xi(1)-ax(1))/(ax(2)-ax(1))*w(3)+w(1),(xi(2)-ax(1))/(ax(2)-ax(1))*w(3)+w(1)];
  yo = [(yi(1)-ax(3))/(ax(4)-ax(3))*w(4)+w(2),(yi(2)-ax(3))/(ax(4)-ax(3))*w(4)+w(2)];
  annotation ("doublearrow",xo,yo,'color',cor);
 endfunction

 function result = clean(x)
    threshold=1e-10;
    result = x;
    if (abs(real(result)) < threshold)  result = 0 + i*imag(result); end
    if (abs(imag(result)) < threshold)  result = real(result); end
endfunction
