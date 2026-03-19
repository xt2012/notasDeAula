prt = false;

function [p2,ordem2,rts]=FatorarIrredutiveis(p_in)
    n=length(p_in);
    n_pol=1;
    while n>3 %repita enquanto ordem do polinomio > 2
          [ps2_out,p_div]=encontrar_divisor_quadratico(p_in);
          [rts(n-1),rts(n-2),delta]=baskara(ps2_out);
          if (delta<0)
              p(n_pol,:)=[1,-2*real(rts(n-1)),abs(rts(n-1))^2];
              ordem(n_pol)=-2;
              n_pol=n_pol+1;
          else
              p(n_pol,:)=[0, 1, -rts(n-1)];
              ordem(n_pol)=1;
              n_pol=n_pol+1
              p(n_pol,:)=[0, 1, -rts(n-2)];
              ordem(n_pol)=1;
              n_pol=n_pol+1;
          end
          p_in=p_div;
          n=n-2;   % diminua em 2 a ordem do polinomio
    end;

    if(n==3)  % ultima equacao ordem 2
        [rts(n-1),rts(n-2)]=baskara(p_in);
        if (delta<0)
              p(n_pol,:)=[1,-2*real(rts(n-1)),abs(rts(n-1))^2];
              ordem(n_pol)=-2;
        else
              p(n_pol,:)=[0, 1, -rts(n-1)];
              ordem(n_pol)=1;
              n_pol=n_pol+1;
              p(n_pol,:)=[0, 1, -rts(n-2)];
              ordem(n_pol)=1;
        end
    else     % ultima equacaoo ordem 1
         rts(n-1)=-p_in(2)/p_in(1);
         p(n_pol,:)=[0, 1, -rts(n-1)];
         ordem(n_pol)=1;
    end;

    for(i=1:n_pol)   % se a fase for pequena, trate como real
       if (ordem(i)==-2)
            r=roots(p(i,:));
            if( atan( abs(imag(r(1))) / abs(real(r(1)))) < 1e-4 )
                p(i,:)=[1, -sign(real(r(1)))*abs(r(1))];
                ordem(i)=1;
                n_pol=n_pol+1;
                p(n_pol,:)==[1, -sign(real(r(1)))*abs(r(1))];
                ordem(n_pol)=1;
            end
        end
     end

    if(n_pol>2)
        for(i=1:n_pol) % raizes reais parecidas sao tratatas com multiplas
            for(j=i+1:n_pol)
                 if ( ordem(i)!=-2 && ordem(j)!=-2)
                     r1=roots(p(i,:));
                     r2=roots(p(j,:));
                     if( abs(r1(1)-r2(1))<1e-4)
                        r = (real(r1(1))+real(r2(1)))/2;
                        p(min(i,j),:)= [0 1 -r];
                        ordem(min(i,j))=ordem(i)+ordem(j)
                        p(max(i,j),:)=1;
                        ordem(max(i,j))=0;
                     end
                 end
            end
        end
    end

    k=1;
    for(i=1:n_pol) % reorganiza o array de saida retirando os eliminados
        if (ordem(i)!=0)
            p2(k,:)=p(i,:);
            ordem2(k)=ordem(i);
            k=k+1;
        end
    end

    [g,k]=sort(real(rts));
    rts = clean(rts(k))';
endfunction

function rts=bairstow_zeros(p_in,prt)
    n=length(p_in);
    while n>3 %repita enquanto ordem do polinomio > 2
          [ps2_out,p_div]=encontrar_divisor_quadratico(p_in);
          [rts(n-1),rts(n-2)]=baskara(ps2_out);
          if (prt)
              disp(p_in);
              disp(ps2_out);
              disp(rts(n-2:n-1)');
              printf("\n");
          end
          p_in=p_div;
          n=n-2;   % diminua em 2 a ordem do polinomio
    end;
    if (prt)  disp(p_in) end
    if(n==3)  % ultima equacao ordem 2
        [rts(n-1),rts(n-2)]=baskara(p_in);
        if (prt) disp(rts(n-2:n-1)') end
    else     % ultima equacaoo ordem 1
         rts(n-1)=-p_in(2)/p_in(1);
         if (prt) disp(rts(n-1)) end
    end;
    [g,k]=sort(real(rts));
    rts = clean(rts(k))';
endfunction

function [u,bq]=encontrar_divisor_quadratico(a)
    u= [1.0,0.1,0.1];
    du=[0.0,0.0,0.0];
    n=length(a);
    for (k=1:500)
       [bq, rb] = fatorar_p2(a,u);
       b=[bq rb];
       [cq, rc] = fatorar_p2(b,u);
       c=[cq rc];
       delta =  c(n-1)*c(n-3)-c(n-2)^2;
       du(2) =  ( b(n)  *c(n-3) - c(n-2)*b(n-1) ) / delta;
       du(3) =  ( c(n-1)*b(n-1) - b(n)  *c(n-2) ) / delta;
       if (norm(du)<1e-14) break end
       u=u+du;
    end
endfunction

function [b_out, r] = fatorar_p2(a, u)
    a=[a 0 0];
    n = length(a);  % fator = s^2 + u(1)*s + u(2)
    b = zeros(1, n);
    for i = 1:n-2
        q = [0 zeros(1, 2)];
        q(1) = a(i) / u(1);
        b(i) = q(1);
        a(i:i+2) = a(i:i+2) - q(1) * u;
    end
    b_out = b(1:n-4);
    r=b(n-3:n-2);
end

function [r1, r2, delta] = baskara(u)
    delta = u(2)^2 - 4 * u(1) * u(3);
    r1 = (-u(2) + sqrt(delta)) / (2 * u(1));
    r2 = (-u(2) - sqrt(delta)) / (2 * u(1));
endfunction

function result = clean(x)
    threshold=1e-10;
    result = x;
    if (abs(real(result)) < threshold)  result = 0 + i*imag(result); end
    if (abs(imag(result)) < threshold)  result = real(result); end
endfunction
