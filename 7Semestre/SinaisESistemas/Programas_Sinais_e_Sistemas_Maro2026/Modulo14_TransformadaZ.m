prt=false

% Degrau unitario
function  x=f_u(t)
   x =double(t>=0);
end;

   % Impulso Discreto
function  x=f_d(n)
   x = double(n==0);
end;

function f_deff=itransformadaZ_Hq(Hq_in, forma2, prt)
    z=poly(0,'z');
    Hz_in=horner(Hq_in,1/z)
    f_deff=itransformadaZ_Hz(Hz_in, forma2)
endfunction

function disp_pfss_Z_fga(Hz)
    [num_c,den_c]=tfdata(Hz,'v');
    [r, p, k, e] = residuez (num_c, den_c);
    z=tf('z',1);
    cmpx=true;
    K=length(k);
    if (K>0)
      D=k(1)*z^(K-1);
      for n=2:K
         if (k(n)!=0) D=D+k(n)*z^(K-n); end
      endfor
      D
    endif
    for n=1:length(r)
      if (abs(imag(p(n)))<1e-8)
            D=tf(real(r(n))/(z-real(p(n)))^e(n))
      else
        cmpx=!cmpx;
        % R1=(A'*(s-r)+A*(s-r'))/((s-r)*(s-r'))
        num=2*real(r(n))* (z-real(p(n)))-2*imag(r(n))*imag(p(n));
        den = (z-real(p(n)))^2+imag(p(n))^2;
        if (cmpx)
          D=tf(num/den)
        endif
      end
    end
endfunction





function f_out=itransformadaZ_Hz(Hz_in,T, forma2)
    z=tf('z',T);
    H=minreal(Hz_in/z);
    [num_c,den_c]=tfdata(H,'v');
    [r, p, k, e] = residue (num_c, den_c);
    cmpx=true;
    fn_d='';
    fn="fn= @(n) (";

    R=length(r);
    K=length(k);
    P=length(p);

   % soma_p = 0;
   % for (n=1:P) soma_p = soma_p + p(n); end
   % soma_p;
   % if  (soma_p==0) */%somente numerador
   if isempty(p)
        fn=Z_inversa_numerador(num_c,P,fn);
    else
        if (K>0) %fracao impropia
            fn=Z_inversa_numerador(k,0,fn);
        end
        for n=1:length(r)
           if (abs(imag(p(n)))<1e-8) %polos reais
              if (e(n)==1)&&(abs(r(n))>1e-8)
                  if (abs(p(n))<1e-5)
                     str=sprintf("%.5f*f_d(n)",r(n));
                   elseif  (abs(p(n)-1)<1e-5)
 					           str=sprintf("%.5f*f_u(n)",r(n));
                   else
                     str=sprintf("%.5f*(%.5f).^n.*f_u(n)",r(n),p(n));
                   end
                   fn=strcat(fn,"+(",str,")");
               elseif (e(n)==2)&&(abs(r(n))>1e-8)
                   if(abs(p(n))<1e-5)
					               str=sprintf("%.5f*f_d(n-1)",r(n));
                    elseif (abs(p(n)-1)<1e-5)
					               str=sprintf("%.5f.*n.*f_u(n)",r(n));
                    else
					               str=sprintf("%.5f*n.*(%.5f).^n.*f_u(n)",r(n)/p(n),p(n));
				  	        end
                    fn=strcat(fn,"+(",str,")");
               else
                   if(abs(r(n))>1e-8)
                      if(abs(p(n))<1e-5)
					                str=sprintf("%.5f*f_d(n-2)\n",r(n));
                      elseif (abs(r(n)-1)<1e-5)
  					               str=sprintf("%.5f.*(n.^2).*f_u(n)",r(n));
                      else
					                 str=sprintf("%.5f*(n.^2).*(%.5f).^n*.f_u(n)",p(n)/2/p(n)^2,p(n));
					            end
                      fn=strcat(fn,"+(",str,")");
                end
              end
            else   % ate aqui polos reais
			    cmpx=!cmpx;
			    if (cmpx)
             if (imag(p(n))> 0)
               p(n) = conj(p(n));
             endif
				     teta = atan2(imag(p(n)),real(p(n)));
				     phi = atan2(imag(r(n)),real(r(n)));
				     P = abs(p(n));
				     M = abs(r(n));
				     if (forma2);
					   if (abs(P-1)<1e-5)
						    str=sprintf("2*(%.5f)* cos(n*(%.5f)-(%.5f)).*f_u(n)\n",M,teta,phi);
					   else
						    str=sprintf("2*(%.5f)*(%.5f).^n .* cos(n*(%.5f)-(%.5f)).*f_u(n)\n",M,P,teta,phi);
					   end
				   else
					 if (abs(P-1)<1e-5)
						 str=sprintf("2*(%.5f)*((%.5f)*cos(n*(%.5f))+(%.5f)*sin(n*(%.5f))).*f_u(n)\n"...
								  ,M,cos(phi),teta,sin(phi),teta);
					  else
						 str=sprintf("2*(%.5f)*(%.5f).^n .* ((%.5f)*cos(n*(%.5f))+(%.5f)*sin(n*(%.5f))).*f_u(n)\n"...
								 ,M,P,cos(phi),teta,sin(phi),teta);
					   end
				   end
				   fn=strcat(fn,"+(",str,")");
				end  % processando segundo polo complexo conjugado
            end  % else polos reais
		end  %para todos os polos
    end % else somente numerador

    fn=strcat(fn,");");
    eval(fn);
    f_out=fn;
endfunction

function f_deff=Z_inversa_numerador(num,idx,f_deff)
        h=num; % intensidade dos impulsos
        N=length(h);
        for (i=1:N )
            n = N-i-idx+1; %adiantamento do inpulso
            if h(i)!=0
                if(n>0)
                   str = sprintf("(%.2f)*f_d(n+%d)",h(i),n);
                elseif(n<0)
                   str = sprintf("(%.2f)*f_d(n%d)",h(i),n);
                else
                   str = sprintf("(%.2f)*f_d(n)",h(i));
                end
                if (abs(h(i))>1e-5)
                  f_deff=strcat(f_deff,"+(",str,")");
                endif
            end
        end
endfunction

%Tentativa de transformar um kernel de convolucao h (resposta ao impulso fir)
% em um filtro recursivo iir somente com polos
function Hq1=DesignMQ(h,K)
    N=length(h)
    A=zeros(N-K,K)
    y=h(K+1:N)'
    A(:,1)=h(K:N-1)'
    for (k=1:K-1)
        A(:,k+1)=h(K-k:N-k-1)'
    end
    a=inv(A'*A)*(A'*y)
    a1 = [1; -a]
    disp(a1)
    disp(h(K-1))
    Hq1 = h(K-1)/poly(a1,'q','coeff')
endfunction


