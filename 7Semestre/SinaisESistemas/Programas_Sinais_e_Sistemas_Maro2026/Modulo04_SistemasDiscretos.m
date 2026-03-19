prt=false;

function plot_vbar(t,x,cor)
  hold on; grid on;
  for i = 1:length(x)
    plot([t(i),t(i)],[0,x(i)],cor,'LineWidth',1)
  endfor
  axis([min(t),max(t)]);
end

function [y_out,n_out]=conv_fga(x,xi,h,prt)
    nx = length(x);
    nh = length(h);
    i_ini = xi - nh +1;
    i_fin = xi + nx + nh - 2;
    ncols = i_fin-i_ini+1;
    nlins = nx + nh -1;;
    indice=[i_ini:i_fin];
    x_ext=zeros(1,ncols);
    for coluna=1:nx
        x_ext(coluna+nh-1)=x(coluna);
    end
    h_r = flip(h);
    h_des=zeros(nlins,ncols);
    for linha=1:nlins
        for i = 1:nh
            h_des(linha,i+linha-1)=h_r(i);
        end
       y_out(linha)=sum(x_ext.*h_des(linha,:));
    end
    for linha=1:nlins
       h_des(linha,ncols+1)=y_out(linha);
       h_des(linha,ncols+2)=xi+linha-1;
    end
    if (prt)
        disp(indice);
        disp(x_ext);
        disp(h_des);
    end
    n_out = h_des(:,ncols+2)';
endfunction

function y_out=rec1A(x_in)
    xi=[0 0 x_in];
    for (i=3:length(xi));
         yi(i)=(xi(i)+xi(i-1)+xi(i-2))/3;
    end
    y_out=yi(3:length(xi));
endfunction

function y_out=rec1B(x_in)
    xi=[0 0 0 x_in];
    for (i=4:length(xi));
         yi(i)=(4*xi(i)+3*xi(i-1)+2*xi(i-2)+xi(i-3))/10;
    end
    y_out=yi(4:length(xi));
endfunction

function y_out=exemplo1(x_in,y_ini)
    xi=[0 0 x_in];
    yi=[ nan  y_ini];
    for (i=3:length(xi))
          yi(i)=5*xi(i)+3*xi(i-2)-3/4*yi(i-1);
    end
    y_out=yi(3:length(xi));
endfunction


function y_out=rec2(x_in,y_ini)
    xi=[0 0 x_in];
    yi=[ nan  y_ini];
    for (i=3:length(xi))
          yi(i)=5*xi(i)+10*xi(i-2)-3/4*yi(i-1);
    end
    y_out=yi(3:length(xi));
endfunction

function y_out=rec3_FIR(x_in)
    xi=[zeros(1,11) x_in];
    for i=12:length(xi)
         yi(i)=mean( xi(i-11:i) );
    end
    y_out=yi(12:length(xi));
endfunction

function y_out=rec3_IIR(x_in,y_ini)
    xi=[zeros(1,12) x_in];
    yi=[zeros(1,11)* nan  y_ini];
    for (i=13:length(xi))
        yi(i)=yi(i-1)+1/12*(xi(i)-xi(i-12));
    end
    y_out=yi(13:length(xi));
endfunction

function y_out=Sistema_LP(x_in,y1,y2,y3)
    xi=[0 0 0  x_in];
    yi=[y3 y2 y1];
    for (i=4:length(xi))
       yi(i)=0.8771*yi(i-1)-0.2455*yi(i-2)+0.02224*yi(i-3)+...
       0.04328*xi(i)+0.1298*xi(i-1)+0.1298*xi(i-2)+0.04328*xi(i-3);
    end
    y_out=yi(4:length(xi));
endfunction

function y_out=Sistema_HP(x_in)
    xi=[zeros(1,5) x_in];
    for (i=6:length(xi))
         yi(i)=   xi(i)  +2*xi(i-1) +4*xi(i-2)...
              -4*xi(i-3)-2*xi(i-4) -xi(i-5);
    end
    y_out=yi(6:length(xi));
endfunction

function y_out=recorrencia_teste(x_in)
    x=[0 0 0 0 x_in]
    for (i=5:length(x))
         y(i)=4*x(i)-x(i-1)+5*x(i-2)+2*x(i-3)+3*x(i-4)
    end
    y_out=y(5:length(x))
endfunction


function y_out=recorrencia2b(x_in,y_ini)
    x=[0 0 x_in]
    y=[ nan  y_ini]
    for (i=3:length(x))
          y(i)=5*x(i)+3*x(i-2)-3/4*y(i-1)
    end
    y_out=y(3:length(x))
endfunction



function y_out=recorrencia3_FIR_centrada(x_in)
    x=[zeros(1:7) x_in]
    for i=7:length(x)-5
         y(i)=mean( x(i-6:i+5) )
    end
    y_out=y(7:length(x)-5)
endfunction


function y_out=Sistema_Eletrico_1aOrdem(x_in,y_ini)
    x=[0 x_in]
    y=[y_ini]
    for (i=2:length(x))
          y(i)=5/6*y(i-1)+1/6*x(i)
    end
    y_out=y(2:length(x))
endfunction

function y_out=Sistema_Eletrico_2aOrdem(x_in,y1,y2)
    x=[0 0 x_in]
    y=[y1 y2]
    for (i=3:length(x))
         y(i)=125/78*y(i-1)-50/78*y(i-2)+3/78*x(i)
    end
    y_out=y(3:length(x))
endfunction

function y_out=Sistema_Mecanico_2aOrdem(x_in,y1,y2)
    x=[0 0 x_in]
    y=[y1 y2]
    for (i=3:length(x))
         y(i)=250/137*y(i-1)-125/137*y(i-2)...
              +12/137*x(i)-10/137*x(i-1)
    end
    y_out=y(3:length(x))
endfunction


function y_out=Prova(x_in,y1,y2)
    x=[0 0 x_in]
    y=[y1 y2]
    for (i=3:length(x))
         y(i)=1.62*y(i-1)-0.67*y(i-2)...
              +0.054*x(i)
    end
    y_out=y(3:length(x))
endfunction

