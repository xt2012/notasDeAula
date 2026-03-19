prt = false;

function [z_m,z_f] = polar_fga(z)
    z_m=abs(z);
    z_f=atan2(imag(z),real(z));
endfunction

function [z_m,z_f] = polar_d_fga(z)
    z_m=abs(z);
    z_f=atan2d(imag(z),real(z));
endfunction
function [zr, zi] = ret_fga(z)
    zr = real(z);
    zi = imag(z);
endfunction

% y=clean(x,1e-14)
function y=cleancomplex(x)
   e=1e-14;
   y=x;
   index=find(abs(x)<e);
   y(index)=0;;
   index_n=find(abs(x)>e);
   index=find(abs(real(x(index_n)))<e);
   y(index)= i *imag(x(index));
   index=find(abs(imag(x(index_n)))<e);;
   y(index)=real(x(index))
endfunction

function yp = polyval_fga(u, xp)
  if(size(u)(1)==1) u=u'; end
  if(size(xp)(1)==1) xp=xp'; end
  N = length(xp);
  M = length(u);
  A = zeros(N, M);
  for i = 1:M
     A(:, M-i+1) = xp.^(i-1);
  end
  yp(1,:)=A*u;
end

function yp = hval_fga(H, xp)
  [num,den] = tfdata(H,'v');
  N=polyval(num,xp);
  D=polyval(den,xp);
  yp=N./D;
end

function Plot_Funcao_Complexa(H,a,b,azi,incl)
    f_mag = @(x,y) abs(hval_fga(H,complex(x,y)));
    f_fase = @(a,b) atan2(imag(hval_fga(H,complex(a,b))),real(hval_fga(H,complex(a,b))));
    [xi,yi]=meshgrid(a,b);
	subplot(121)
  	zi = f_mag(xi,yi);
    hold on; grid on;
	  mesh(xi,yi,zi);
	  colormap(jet);
	  xlabel('real(s)');
    ylabel('imag(s)');
    zlabel('modulo(z)');
    view(azi,incl);
   subplot(122)
    zi = f_fase(xi,yi);
    hold on; grid on;
	  mesh(xi,yi,zi);
	  colormap(jet);
	  xlabel('real(s)');
    ylabel('imag(s)');
    zlabel('fase(z)');
    view(azi,incl);
end

function [R,P]=residuos(H)
  [num_c,den_c]=tfdata(H,'v');
   [r, p, k, e] = residue (num_c, den_c);
   index = find(e==1);
   R=r(index);
   P=p(index);
endfunction


function disp_pfss_fga(H)
    [num_c,den_c]=tfdata(H,'v');
    [r, p, k, e] = residue (num_c, den_c);
    s=tf('s');
    cmpx=true;
    K=length(k);
    if (K>0)
      D=k(1)*s^(K-1);
      for n=2:K
         if (k(n)!=0) D=D+k(n)*s^(K-n); end
      endfor
      D
    endif

    for n=1:length(r)
      if (abs(imag(p(n)))<1e-8)
            D=tf(real(r(n))/(s-real(p(n)))^e(n))
      else
        cmpx=!cmpx;
        % R1=(A'*(s-r)+A*(s-r'))/((s-r)*(s-r'))
        num=2*real(r(n))* (s-real(p(n)))-2*imag(r(n))*imag(p(n));
        den = (s-real(p(n)))^2+imag(p(n))^2;
        if (cmpx) D=tf(num/den) end
      end
    end
endfunction

function H=Funcao_de_Transferencia(a,b)
  npolos=length(a)-1
  nzeros=length(b)-1
  q=poly(0,'q')
  num=0;
  den=1;
  for p=1:npolos
    den=den+a(p+1)*q^p;
  end
  for p=0:nzeros
     num=num+b(p+1)*q^p;
  end
  H=num/den;
endfunction

function plot_polos_zeros(Hs)
   pol=pole(Hs);
   zer=zero(Hs);
   hold on; grid on;
   title('Posicao dos polos e zeros');
   plot(real(pol),imag(pol),'x');
   plot(real(zer),imag(zer),'o');
   angle = linspace(0, 2* pi , 1024);
   r=max(abs(pol));
   plot(r*cos(angle),r*sin(angle),'r--');
   lim=1.05*max(abs([pol;zer]));
   axis ([-lim,lim,-lim,lim]);
endfunction

