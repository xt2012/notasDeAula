prt=false;

function y_out=dx(t,x)
    N=length(t);
    y_out(1)=(x(t(2))-x(t(1)))/(t(2)-t(1));
    for(n=2:N-1)
       y_out(n)=(x(t(n+1))-x(t(n-1)))/(t(n+1)-t(n-1));
    end
    y_out(N)=(x(t(N))-x(t(N-1)))/(t(N)-t(N-1));
endfunction

function y_out=intx(t,x)
   N=length(t);
   y_out(1)=0;
   for(n=2:N)
       y_out(n)=y_out(n-1) + (x(t(n))+x(t(n-1)))/2 * (t(n)-t(n-1));
   end
endfunction

function y_out=rec1(t,x,y_ini)
    N=length(t);
    y_out=[y_ini];
    for (i=2:N)
          y_out(i)=0.62*x(t(i))+0.38*y_out(i-1);
    end
endfunction

function y_out=rec1b(t,x,y_ini)
    N=length(t);
    y_out=[y_ini]
    for (i=2:N)
          y_out(i)=0.38*x(t(i))+0.62*y(i-1)
    end
endfunction


