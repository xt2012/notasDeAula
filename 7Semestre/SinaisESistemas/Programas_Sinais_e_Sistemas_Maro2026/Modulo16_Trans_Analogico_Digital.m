prt=false;

function Hz=Hs2Hz_elementos(Hs,T)
    [s_num,s_den] = tfdata(Hs,'v');
    if (length(s_num)<5) s_num =[zeros(1,5-length(s_num)),s_num]; end
    if (length(s_den)<5) s_den =[zeros(1,5-length(s_den)),s_den]; end
    q=tf('q',T);
    G1= tf(1,q);
    num=0;
    if (s_num(5)!=0)  num += G1*s_num(5); end
    if (s_num(4)!=0)  num+=(G1-q)/T*s_num(4); end
    if (s_num(3)!=0)  num+=(G1-2*q+q^2)/T^2*s_num(3); end
    if (s_num(2)!=0)  num+=(G1-3*q+3*q^2-q^3)/T^3*s_num(2); end
    if (s_num(1)!=0)  num+=(G1-4*q+6*q^2-4*q^3+q^4)/T^4*s_num(1); end
    den=0;
    if (s_den(5)!=0)  den+= G1*s_den(5); end
    if (s_den(4)!=0)  den+=(G1-q)/T*s_den(4); end
    if (s_den(3)!=0)  den+=(G1-2*q+q^2)/T^2*s_den(3); end
    if (s_den(2)!=0)  den+=(G1-3*q+3*q^2-q^3)/T^3*s_den(2); end
    if (s_den(1)!=0)  den+=(G1-4*q+6*q^2-4*q^3+q^4)/T^4*s_den(1); end
    Hq=(num/den);
    [s_num,s_den] = tfdata(Hq,'v');
    Hq=tf((s_num/s_den(end)),(s_den/s_den(end)),T,"Variable","q");
    Hz=Hq2Hz(Hq,2*pi/T);
endfunction

function Hz=HsHz_complex(Hs,T) %par complexo
    z=tf('z',T);
    [zer,pol] = tfdata(Hs,'v');
    if length(zer)==1    zer=[0 zer(1)]; end
    B = zer(end-1);
    C = zer(end);
    a =  pol(2)/(2*pol(1));
    wd = sqrt(pol(3)/pol(1) - a^2);
    num = B*z*(z-exp(-a*T)*cos(wd*T))+...
     (C-a*B)/wd*z*exp(-a*T)*sin(wd*T);
    den = z^2-2*z*exp(-a*T)*cos(wd*T)+exp(-2*a*T);
    Hz=T*num/den;
endfunction

function Hz=HsHz_simples(Hs,T)
    z=tf('z',T);%polo simples
    [zer,pol] = tfdata(Hs,'v');
    A = zer(end);   % sempre o coeficiente correto
    a = pol(end)/pol(1);
    num = A*z;
    den = z-exp(-a*T);
    Hz=T*num/den;
endfunction

function Hz=HsHz_duplo(Hs1,Hs2,T)
   z=tf('z',T); % polo duplo
   [zer,pol] = tfdata(Hs1,'v');
    A1 = zer(end);   % sempre o coeficiente correto
    a = pol(end)/pol(1);
   [zer,pol] = tfdata(Hs2,'v');
   A2 = zer(end);
   den1 = z-exp(-a*T);
   num1 = A1*z;
   num2 = A2*z*T*exp(-a*T);
   num  = num1*(den1) + num2;
   Hz=T*num/den1^2;
endfunction

function Hz=HsHz_triplo(Hs1,Hs2,Hs3,T)
   z=tf('z',T);%polo triplo
   [zer,pol] = tfdata(Hs1,'v');;
   A1 = zer(end);
   a = pol(end)/pol(1);
   [zer,pol] = tfdata(Hs2,'v');;
   A2 = zer(end);
   [zer,pol] = tfdata(Hs3,'v');;
   A3 = zer(end);
   den1 = z-exp(-a*T);
   num1 = A1*z;
   num2 = A2*z*T*exp(-a*T);
   num3 = 1/2*A3*z^2*T^2*exp(-a*T)*(z+exp(-a*T));
   num  = num1*den1^2 + num2 * den1 + num3;
   Hz=T*num/den1^3;
endfunction


function  Ht = soma_tf_hz(H1,H2,T)
  warning('off', 'all');
  [num1,den1]=tfdata(H1,'v');
  [num2,den2]=tfdata(H2,'v');

  syms zz;
	Hzz1 = poly2sym(num1,'zz')/poly2sym(den1,'zz');
	Hzz2 = poly2sym(num2,'zz')/poly2sym(den2,'zz');
	Hzz3=simplify(Hzz1+Hzz2);
	[n,d]=numden(Hzz3);

  num3=sym2poly(n);
  den3=sym2poly(d);
	Ht = tf(num3/den3(1),den3/den3(1),T);
endfunction

function Hz=Hs2Hz_invariancia(H,T)
    [num_c,den_c]=tfdata(H,'v');
    [r, p, k, e] = residue (num_c, den_c);
    s=tf('s');
    cmpx=true;
    Hz = tf(0);
    if (length(k)!=0)
      if (k(end)!=0)
          Hz=k(end)*z;
      end
    end
    for n=length(r):-1:1
      if (abs(imag(p(n)))<1e-8)
            if (e(n)==1 && n==length(r))
               D=tf(real(r(n))/(s-real(p(n))));
               Hza=HsHz_simples(D,T);
               Hz = soma_tf_hz(Hz,Hza,T);
            elseif (e(n)==1&& e(n+1)==1)
               D=tf(real(r(n))/(s-real(p(n))));
               Hza=HsHz_simples(D,T);
               Hz = soma_tf_hz(Hz,Hza,T);
            elseif (e(n)==2)
               D1=tf(real(r(n-1))/(s-real(p(n))));
               D2=tf(real(r(n  ))/(s-real(p(n)))^2);
               Hza=HsHz_duplo(D1,D2,T);
               Hz = soma_tf_hz(Hz,Hza,T);
            elseif (e(n)==3)
               D1=tf(real(r(n-2))/(s-real(p(n)))  );
               D2=tf(real(r(n-1))/(s-real(p(n)))^2);
               D3=tf(real(r(n  ))/(s-real(p(n)))^3);
               Hza=HsHz_triplo(D,T);
               Hz = soma_tf_hz(Hz,Hza,T);
            end
      else
        cmpx=!cmpx;
        num=2*real(r(n))* (s-real(p(n)))-2*imag(r(n))*imag(p(n));
        den = (s-real(p(n)))^2+imag(p(n))^2;
        if (cmpx)
          D=tf(num/den);
          Hza=HsHz_complex(D,T);
          Hz = soma_tf_hz(Hz,Hza,T);
        end
      end
     end
endfunction

function [ tfobj ] = sym2tf( symobj, Ts)
  % SYM2TF convert symbolic math rationals to transfer function
  if isnumeric(symobj)
      tfobj=symobj;
      return;
  end
  [n,d]=numden(symobj);
  num=sym2poly(n);
  den=sym2poly(d);
  if nargin==1
      tfobj=tf(num,den);
  else
      tfobj=tf(num,den,Ts);
  end
endfunction

function [Hss]= tf2sym( Hs)
  [Num,Den] = tfdata(Hs,'v');
   syms ss
   Hss=poly2sym(Num,ss)/poly2sym(Den,ss)
endfunction

function [funcobj]= sym2func( symobj)
  funcobj = matlabFunction(symobj);
endfunction

function Hz= Hs2Hz_bilinear(Hs,w0,ws)
    %de tf (Hs) para syms (Hss)
    [Num,Den] = tfdata(Hs,'v');
    syms ss;
    Hss=poly2sym(floor(Num),ss)/poly2sym(floor(Den),ss);
    %pre-warp
    wa=ws/pi*tan(pi*w0/ws);
    Hssa = subs(Hss,ss,ss*w0/wa);
    %Transformacao Bilinear
    syms zz;
    Hzz = simplify( subs(Hssa,ss,ws/pi*(zz-1)/(zz+1)) );
    %de syms (Hzz) para tf (Hz)
    [n,d]=numden(Hzz);
    num=sym2poly(n);
    den=sym2poly(d);
    Hz=tf(num/den(1),den/den(1),2*pi/ws);
endfunction


function Hq= Hs2Hq_bilinear(Hs,w0,ws)
  %transformar de tf (Hs)  para tf (Hz)
  Hz= Hs2Hz_bilinear(Hs,w0,ws);
  Hq=Hz2Hq(Hz,ws);
endfunction

function Hq=Hz2Hq(Hz,ws)
  %tranformar de tf (Hz) para syms (Hzz)
  [Num,Den] = tfdata(Hz,'v');
  syms zz;
  Hzz=poly2sym(Num,zz)/poly2sym(Den,zz);

  %horner de zz para 1/qq
  syms qq;
  Hqq = subs(Hzz,zz,1/qq);
  Hqq_s=simplify(Hqq);

  %Transformar de sys (Hqq_s) para tf (Hq)
  [n,d]=numden(Hqq_s);
  num=sym2poly(n);
  den=sym2poly(d);
  Hq=tf(num/den(end),den/den(end),2*pi/ws,"Variable","q");
endfunction

function Hz=Hq2Hz(Hq,ws)
  warning('off', 'all');
  %tranformar de tf (Hz) para syms (Hzz)
  [Num,Den] = tfdata(Hq,'v');
  syms qq;
  Hqq=poly2sym(Num,qq)/poly2sym(Den,qq);

  %horner de zz para 1/qq
  syms zz;
  Hzz = subs(Hqq,qq,1/zz);
  Hzz_s=simplify(Hzz);

  %Transformar de sys (Hqq_s) para tf (Hq)
  [n,d]=numden(Hzz_s);
  num=sym2poly(n);
  den=sym2poly(d);
  z=tf('z',2*pi/ws);
  Hz=tf(num/den(1),den/den(1),2*pi/ws);
endfunction

function Resposta_Frequencia_HsHz(Hs,Hz,w_max,ws_plot)
    if (w_max==0)  w_max=10*max(abs(roots(Hs.den))); end
    w=linspace(-w_max,w_max,1024);
    Hw1=hval_fga(Hs,i*w);
    plot(w,abs(Hw1),'r--');
    hold on; grid on;

   w=linspace(-pi , pi ,1024);
   Hw=hval_fga(Hz,exp( i *w));
   plot(w/(2* pi )*ws_plot,abs(Hw),'b')
 %  axis ([min(w/(2* pi )*ws_plot),max(w/(2* pi )*ws_plot),0,max(abs(Hw))]);;
   axis([-w_max,w_max,0,max([abs(Hw),abs(Hw1)])]);
   xticks (linspace(-w_max,w_max,9));
   ax1=gca;
   set(ax1, 'XAxisLocation', 'top','YAxisLocation','Right');
   ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
   xticks (ax2,[-pi:pi/4:pi]/(2*pi)*ws_plot);
   set(ax2, 'XAxisLocation', 'bottom','YAxisLocation','Left');
   set(ax2, 'XLim', get(ax1, 'XLim'));
   set(ax2, 'XTickLabel', {'-pi','-3pi/4','-pi/2','-pi/4','0','pi/4','pi/2','3pi/4','pi'});
   grid on;


endfunction

