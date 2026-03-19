prt=false

function y_out=Sistema_Eletrico_1aOrdem(x_in,y_ini)
    x=[0 x_in];
    y=[y_ini];
    for (i=2:length(x))
          y(i)=5/6*y(i-1)+1/6*x(i);
    end
    y_out=y(2:length(x));
endfunction

function y2=applyIIR(Hz,x,xp,yp)
   [b,a]=tfdata(Hz,'v');
   npolos=length(a)-1;
   nzeros=length(b)-1;
   N=length(x);
   for n=1:N
     y2(n)=0;
     for np=1:npolos
        index=n-np;
        if index<=0
           y2(n)=y2(n)-a(np+1)*yp(npolos+index);
        else
           y2(n)=y2(n)-a(np+1)*y2(index);
        end
     end
     for nz=0:nzeros
        index=n-nz;
        if index<=0
            y2(n)=y2(n)+b(nz+1)*xp(nzeros+index);
        else
            y2(n)=y2(n)+b(nz+1)*x(n-nz);
        end
     end
    end
endfunction

function Resposta_Frequencia_Hz(Hz,ws_plot)
 w=linspace(- pi , pi ,1024);
 Hw=hval_fga(Hz,exp( i *w));
  subplot(212)
   plot(w/(2* pi )*ws_plot,atan2(imag(Hw),real(Hw)));
   axis ([min(w/(2* pi )*ws_plot),max(w/(2* pi )*ws_plot),-pi,pi]);;
   title("Fase em Radianos")
   grid on
 subplot(211)
   plot(w/(2* pi )*ws_plot,abs(Hw))
   axis ([min(w/(2* pi )*ws_plot),max(w/(2* pi )*ws_plot),0,max(abs(Hw))]);;
   title("|H(jw)|")
   ax1=gca;
   set(ax1, 'XAxisLocation', 'top','YAxisLocation','Right');
   ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
   xticks (ax2,[-pi:pi/4:pi]/(2*pi)*ws_plot);
   set(ax2, 'XAxisLocation', 'bottom','YAxisLocation','Left');
   set(ax2, 'XLim', get(ax1, 'XLim'));
   set(ax2, 'XTickLabel', {'-pi','-3pi/4','-pi/2','-pi/4','0','pi/4','pi/2','3pi/4','pi'});
   grid on;

endfunction


function DiagramaBode_Hz(Hz,logf1,logf2,cor)
   logf=linspace(logf1,logf2,1024);
   w=2* pi *(10.^logf);
   Hw=horner(Hq,exp( i *w))
  subplot(211)
   plot(logf,20*log10(abs(Hw)),cor)
   xtitle("Amplitude em dB","log10(f) Hertz","|H(jw)|dB")
  subplot(212)
   plot(logf,atan(imag(Hw),real(Hw)),cor);
   xtitle("Fase","log10(f) Hertz","Fase(rad)")
endfunction

function y=SomaResiduos(Hz_in,n)
  Hz=Hz_in
  den=Hz.den;
  num=Hz.num;
  [ps,ordem]=FatoraTermosIrredutiveis(horner(Hz.den,s)) % troca de z para s
  pz=horner(ps,z) % volta de s para z
  ordem=coeff(ordem) % 1-simples 2-duplo 3-triplo -2-complexconj
  N=length(pz)
  y=n*0
  for(i=1:N)
      den_res=1
      for(j=1:N) % denominador para o caxis (lculo dos residuos
          if(i!=j)
               if(ordem(j)>1)      den_res=den_res*pz(j)^ordem(j)
               else                   den_res=den_res*pz(j)   end
          end
      end
      if (ordem(i)==1)    %1 polo simples
          HA = num/den_res * z^(n-1)
          polo=-(coeff(pz(i))(1))
          A=horner(HA,polo)
          y=y+A
      elseif    (ordem(i)==2) %1 polo duplo
          HA = num/den_res * z^(n-1)
          polo=-(coeff(pz(i))(1))
          HAd=derivat(HA)
          A=horner(HAd,polo)
          y=y+A
      elseif    (ordem(i)==3) % 1 polo triplo
          HA = num/den_res * z^(n-1)
          polo=-(coeff(pz(i))(1))
          HAd=derivat(HA)
          HAdd=derivat(HAd)
          A=1/2*horner(HAdd,polo)
          y=y+A
      elseif    (ordem(i)==-2) % 2 polos complexos conjugados
          polo=baskara(horner(pz(i),s)) %troca z por s
          HA = num/(den_res*(z-conj(polo))) * z^(n-1)
          A=horner(HA,polo)
          y=y+2*real(A)
      end
  end
endfunction

