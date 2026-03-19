prt=false
pkg load signal
pkg load control
pkg load symbolic
s=tf('s');

% Impulso Discreto
function  x=f_d(n)
   x = double(n==0);
end;

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

% Funcao Quadratica
function  x=f_q(t)
   x =t.^2.*f_u(t);
end;

function ft=ilaplace_fga(H,tipo1)

    s=tf('s');
    if isnumeric(H) H=tf(1); end
    H=minreal(H);
    [num_c,den_c]=tfdata(H,'v');

    [r, p, k, e] = residue (num_c, den_c);
    cmpx=true;
    f_deff_d='';
    f_deff="ft=@(t)(";

    if (!isempty(k)) % fracao impropia
        f_deff_d = "f_total(t)=";
        k=flip(k);
        for (i=1:length(k))
            if k(i)!=0
                 strd=sprintf("(%.2f)*f_d",k(i));
                 if i>1
                     for j=1:i-1
                         strd=strcat(strd ,"'");
                      end
                 end
                 strd=strcat(strd , "(t)");
                 if (i>1)
                   f_deff_d = strcat(f_deff_d,"+",strd);
                 else
                   f_deff_d = strcat(f_deff_d,strd);
                 end
            end
        end
    end

    for n=1:length(r)
      if (abs(imag(p(n)))<1e-8)
        if (e(n)==1)
              if(abs(r(n))>1e-8)
                 if(abs(p(n))<1e-5)
                      str=sprintf("%.5f",r(n));
                 else
                      str=sprintf("%.5f*exp(%.5f*t)",r(n),p(n));
                 end
                 f_deff=strcat(f_deff,"+(",str,")");
              end
         elseif (e(n)==2)
              if(abs(r(n))>1e-5)
                 if(abs(p(n))<1e-5)
                       str=sprintf("%.5f*t",r(n));
                  else
                       str=sprintf("%.5f*t.*exp(%.5f*t)",r(n),p(n));
                  end
                 f_deff=strcat(f_deff,"+(",str,")");
              end
          elseif (e(n)==3)
             if(abs(r(n))>1e-8)
                  if(abs(p(n))<1e-5)
                      str=sprintf("%.5f*1/2*t.^2",r(n));
                  else
                      str=sprintf("%.5f*1/2*t.^2.*exp(%.5f*t)",r(n),p(n));
                  end
                  f_deff=strcat(f_deff,"+(",str,")");
              end
          end
      else
        cmpx=!cmpx;
        if (cmpx)
              alpha = real(p(n));
              w = imag(p(n)');
              a = real(r(n));
              if imag(p(n))< 0
                   b = -imag(r(n));
              else
                   b = imag(r(n));
              end
              atan2(b,a);
              if(abs(alpha)<1e-4)
                   str=sprintf("2*(%.5f*cos(%.5f*t)-(%.5f)*sin(%.5f*t))",a,w,b,w);
                   str2=sprintf("2*(%.5f)*cos(%.5f*t+(%.5f))",sqrt(a^2+b^2),w,atan2(b,a));
              else
                   str=sprintf("2*exp(%.5f*t).*[%.5f*cos(%.5f*t)-(%.5f)*sin(%.5f*t)]"
                   ,alpha,a,w,b,w);
                   str2=sprintf("2*(%.5f).*exp(%.5f*t).*cos(%.5f*t+(%.5f))",sqrt(a^2+b^2)
                   ,alpha,w,atan2(b,a));
              end
              if tipo1
                  f_deff=strcat(f_deff,"+(",str2,")");
              else
                  f_deff=strcat(f_deff,"+(",str,")");
              end;
        endif
      end
    end

    if (length(f_deff)<=8)
      disp(f_deff_d)
      ft= @(t)  double(n==0);
    else
      if (length(f_deff_d)>0)
        f_deff_d=strcat(f_deff_d,"+ft(t)");
        disp(f_deff_d)
      end
      f_deff=strcat(f_deff," ).*f_u(t);");
       format rat
      disp(f_deff)
      format rat;
      eval(f_deff)
    end
endfunction

function [yo,t,ft]=Plot_Transformada_Inversa(Hs,L,cor,tipo1)
    s=tf('s');
    if isnumeric(Hs) Hs=tf(1); end
    Hs=minreal(Hs);
    [num_c,den_c]=tfdata(Hs,'v');
    [r, p, k, e] = residue (num_c, den_c);
    ft=ilaplace_fga(Hs,tipo1);
    index=find(abs(p)>0);
    min_p=0;
    if( !isempty(index) )
		min_p = min(abs(p(index)));
    end
    if(min_p>0)
       ct = 1/min_p; %constante de tempo
       max_t = 10 * L * ct;
    else
       ct=0.1;
       max_t = 10000*ct;
    end
    plot_impulse=false;
    t=linspace(-ct,max_t,10001);;
    hold on; grid on;
    if (num_c==[1] && den_c==[1])
      clf;
       yo=f_d(t);
       plot(t,yo,cor);
       plot_arrow_fga([0,0],[0,1],'r')
    else
      if isequal(axis, [-0.2000   1.0000  -1.0000   1.0000] )
         clf;
         plot_impulse=true;
      end;
      p=abs(p);
      ws = 50*max(abs(p));
      T=2* pi /ws;
      min_p=0;
      ft;
      yo=ft(t);
      plot(t,yo,cor);
      title("y(t) - Transformada inversa de Laplace de Y(s)");
      hold on; grid on;
      if (!isempty(k)) % fracao impropia
          k=flip(k);
          I=max(abs(k));
          if length(k)>=1
              plot([0,0],[0,max(yo)*k(1)/I]);
              plot_arrow_fga([0,0],[0,max(yo)*k(1)/I],cor);
          end
          if length(k)>=2
              dt=max(t)-min(t);
              plot([-dt/100,-dt/100],[-max(yo)/2*k(2)/I,max(yo)/2*k(2)/I]);
              plot_darrow_fga([-dt/100,-dt/100],[-max(yo)/2*k(2)/I,max(yo)/2*k(2)/I],cor);
          end
          if (plot_impulse)
              plot_arrow_fga([0,0],[0,max(yo)/I],'r')
          end
      else
          if (plot_impulse)
              plot_arrow_fga([0,0],[0,3/4*max(yo)],'r')
          end
      end
    end
endfunction

function  Xs=EntradaNula(Hs,y_ini)
    [num_c,den_c]=tfdata(Hs,'v');
    num_c=flip(num_c);
    den_c = flip(den_c);
    ny = length(y_ini);
    ncond=length(den_c)-1;
    y=y_ini;
    if ny < ncond
        y(1:ncond)=0;
        y(1:ny)=y_ini(1:ny);
    end
    ncond=length(den_c)-1;
    if length(y) < length(den_c)-1
        printf("Entre com %d condicoes iniciais",ncond);
        return
    end
    p = den_c;
    Xs=0;
    s=tf('s');
    for k=1:length(y)
        for coef=k+1:length(p)
            Xs=Xs+y(k)*p(coef)*s^(coef-k-1)
        end
    end

endfunction

function  Ht = soma_tf(H1,H2)
  warning('off', 'all');
  [num1,den1]=tfdata(H1,'v');
  [num2,den2]=tfdata(H2,'v');

  syms ss;
	Hs1 = poly2sym(num1,'ss')/poly2sym(den1,'ss');
	Hs2 = poly2sym(num2,'ss')/poly2sym(den2,'ss');
	Hs3=simplify(Hs1+Hs2);
	[n,d]=numden(Hs3);

  num3=sym2poly(n);
  den3=sym2poly(d);
	Ht = tf(num3,den3);
endfunction

function RespostaSistemaHs(Hs,Xs,y_ini,L,tipo1)
    clf()
    s=tf('s');
    if isnumeric(Hs) Hs=tf(Hs); end
    if isnumeric(Xs) Xs=tf(Xs); end
    hold on;
    plot([0,1e-10],[0,0],'r');
    plot([0,1e-10],[0,0],'b');
    grid on;
    [num_c,den_c]=tfdata(Hs,'v');
    Ys = Hs* Xs;
    Xn= EntradaNula(Hs,y_ini);
    Hn = tf(1,den_c);
    Yn = Hn * Xn;
    if !all(y_ini==0)
      Yt = soma_tf(Ys,Yn);
    else
      Yt=Ys;
    endif
    printf("%s","--Saida--\n"); % linha azul
    [y,t,ft]=Plot_Transformada_Inversa(Yt, L,'b',tipo1);
    printf("\n%s","--Entrada--") ;%linha vermelho pontilhada
    [numx_c,denx_c]=tfdata(Xs,'v');
    I=max(numx_c);
    dt=max(t)-min(t);
    dy=max(y)-min(y);
    if (length(denx_c)==1) %entrada impulso ou double
        if(length(numx_c)==1)
          hold on; grid on;
          if (numx_c(1)!=0)
            plot_arrow_fga([0,0],[0,3/4*max(y)],'r')
          else
             plot(t,t*0,'r');
          endif
          printf("\nfx=@(t)(y=(%.5f)*f_d(t));\n",numx_c(1));
        else
          hold on; grid on;
          if (numx_c(1)!=0)
            plot_arrow_fga([0,0],[0,3/4*max(y)],'r')
          else
            if (numx_c(2)!=0)
                plot_darrow_fga([0,0],[0,3/4*max(y)],'r')
            else
                plot(t,t*0,'r');
            end
          endif
          printf("\nfx=@(t)(y=(%.5f)*f_d(t)+(%.5f)*f_d(t)');\n",numx_c(1),numx_c(2));
        end

 %   elseif (Xs.den==s) % entrada degrau + (impulso)
       elseif isequal(denx_c, [1 0] )
        if(length(numx_c)==1)
          A=numx_c(1);
          B=0;
        else
          A=numx_c(2);
          B=numx_c(1);
          plot_arrow_fga([0,0],[0,3/4*max(y)],'r')
        end
        plot([-max(t)/10,0,0],[0,0,A],"r--");
        plot(t,A*f_u(t),"r--");
        printf("\nfx=@(t)(y=(%.5f)*f_u(t)+(%.5f)*f_d(t));\n",A,B);
 %   elseif (Xs.den==s^2) % entrada rampa (+degrau+impulso)
        elseif isequal(denx_c, [1 0 0] )
        A=numx_c(1);
        B=0;
        C=0;
        if(length(numx_c)==1)
            plot(t,A*f_r(t),"r--");
        elseif(length(numx_c)>=2)
            B=numx_c(2);;
            plot(t,A*f_r(t)+B*f_u(t),"r--");
        end
        if(length(numx_c)==3)
            C=numx_c(3);;
            plot_arrow_fga([0,0],[0,3/4*max(y)],'r')
        end
        printf("\nfx=@(t)(y=((%.5f)*t+(%.5f)).*f_u(t)+(%.5f)*f_d(t));\n",A,B,C)
  %  elseif (Xs.den==s^3) % entrada quadratica (+rampa+degrau+impulso)
        elseif isequal(denx_c, [1 0 0 0] )
        A=numx_c(1);
        B=0;
        C=0;
        D=0;
        if(length(numx_c)==1)
             plot(t,A*f_q(t),"r--");
        elseif(length(numx_c)>=2)
             B=coeff(Xs.num)(2);;
             plot(t,A*f_q(t)+B*f_r(t),"r--");
        elseif(length(numx_c)>=3)
             C=coeff(Xs.num)(3);;
             plot(t,A*f_q(t)+B*f_r(t)+C*f_u(t),"r--");
         end
         if(length(numx_c)==4);
             D=coeff(Xs.num)(4);;
             plot_arrow_fga([0,0],[0,3/4*max(y)],'r');
         end
        printf("\nfx=@(t)(y=((%.5f).*t.^2+(%.5f)*t+(%.5f)).*f_u(t)+(%.5f)*d(t))\n",A,B,C,D)
    else  %demais entradas
         printf("\n");
      %   fx=ilaplace_fga(Xs,tipo1);
         [x_in2,t_in2,fx]=Plot_Transformada_Inversa(Xs,L,'r--',tipo1);
         t_extra=linspace(max(t_in2),max(t),10000);
         plot(t_extra,fx(t_extra),'r--')
    %     plot(t,y,'b')
  end
  legend("entrada x(t)","saida y(t)");
  %  xgrid()
   % gca().data_bounds=[-max(t)/10,min(y)-max(y)/10;max(t),1.1*max(y)];
end


function formula = get_inverse_laplace_formula(H)
  % This function calculates the inverse Laplace transform formula
  % for a given transfer function using symbolic math
  %
  % Input:
  %   H - Transfer function object (created with tf command)
  %
  % Output:
  %   formula - String containing the analytical expression

  % Check and install symbolic package if needed
  try
    pkg list symbolic > /dev/null;
  catch
    % Symbolic package not found, try to install it
    disp('Symbolic package not found. Attempting to install...');
    try
      pkg install -forge symbolic
      pkg load symbolic
      disp('Symbolic package installed successfully.');
    catch installation_error
      warning('Could not install symbolic package. Using residue method instead.');
      % Fall back to residue method
      return_residue_formula = true;
    end
  end

  % Try to load symbolic if it was already installed
  try
    pkg load symbolic
    return_residue_formula = false;
  catch
    return_residue_formula = true;
  end

  % If we couldn't load symbolic, use residue method
  if return_residue_formula
    formula = get_residue_formula(H);
    return;
  end

  % Get transfer function data
  [num, den] = tfdata(H, 'v');

  try
    % Convert to symbolic
    syms s t;

    % Convert numerator and denominator to symbolic polynomials
    num_sym = poly2sym(num, s);
    den_sym = poly2sym(den, s);

    % Create symbolic transfer function
    H_sym = num_sym / den_sym;

    % Calculate inverse Laplace transform
    y_t = ilaplace(H_sym, s, t);

    % Convert to string for display
    formula_raw = char(y_t);

    % Clean up the formula for display
    formula_raw = strrep(formula_raw, '*', ' ');

    % Add Heaviside step function for causal systems if not already there
    if ~contains(formula_raw, 'heaviside')
      formula = ['y(t) = ', formula_raw, ' ?? u(t)'];
    else
      formula = ['y(t) = ', formula_raw];
    end

  catch exception
    % If symbolic approach fails, fall back to residue method
    formula = get_residue_formula(H);
  end
end

function formula = get_residue_formula(H)
  % Fallback function using residue method

  % Get transfer function data
  [num, den] = tfdata(H, 'v');

  try
    % Calculate partial fraction expansion
    [r, p, k] = residue(num, den);

    % Format expression for impulse response
    expr = '';
    for i = 1:length(r)
      if i > 1
        expr = [expr ' + '];
      endif

      if isreal(r(i)) && isreal(p(i))
        expr = [expr, num2str(r(i), '%.4g'), ' e^{', num2str(p(i), '%.4g'), 't}'];
      else
        % Handle complex numbers
        if isreal(r(i))
          r_str = num2str(r(i), '%.4g');
        else
          r_str = ['(', num2str(real(r(i)), '%.4g'), '+', num2str(imag(r(i)), '%.4g'), 'j)'];
        endif

        if isreal(p(i))
          p_str = num2str(p(i), '%.4g');
        else
          p_str = ['(', num2str(real(p(i)), '%.4g'), '+', num2str(imag(p(i)), '%.4g'), 'j)'];
        endif

        expr = [expr, r_str, ' e^{', p_str, 't}'];
      endif
    endfor

    % Handle direct term (k)
    if ~isempty(k) && any(k ~= 0)
      if ~isempty(expr)
        expr = [expr, ' + '];
      endif

      k_terms = '';
      for i = 1:length(k)
        if k(i) ~= 0
          if i > 1
            k_terms = [k_terms, ' + '];
          endif

          if i == 1
            k_terms = [k_terms, num2str(k(i), '%.4g'), ' ??(t)'];
          else
            k_terms = [k_terms, num2str(k(i), '%.4g'), ' ??^(', num2str(i-1), ')(t)'];
          endif
        endif
      endfor

      expr = [expr, k_terms];
    endif

    formula = ['y(t) = ', expr, ' ?? u(t)'];
  catch
    formula = 'Could not calculate analytical form';
  end
end


function ft=ilaplace_fga_s(H)
    s=tf('s');
    if isnumeric(H) H=tf(1); end
    H=minreal(H);
    [num_c,den_c]=tfdata(H,'v');
    Hs1 = poly2sym(num_c,'s')/poly2sym(den_c,'s');
    syms ts;
    fs=ilaplace(Hs1);
    ft = function_handle(fs);
end

function F=laplace_fga_s(f)
   pkg load symbolic;
   pkg load control
   syms ts;
   fs = symfun(f(ts), ts);
   Fs = laplace(fs);
   [ns,ds]=numden(Fs);
   num=sym2poly(ns);
   den=sym2poly(ds);
	 F = tf(num,den);
end
