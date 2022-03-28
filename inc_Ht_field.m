
function [Hx,Hy] = inc_Ht_field(source,test,n,k0,c,OMEGA,direction)

% xsource = source(1,:);
% ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

shift_fact_x = 1i*k0*OMEGA/c.*xtest;
shift_fact_y = 1i*k0*OMEGA/c.*ytest;


switch direction
    
    case 1
        dx = 1i*k0*n*exp(1i*k0*n*ytest);
        dy = 0;
    case -1    
        dx = -1i*k0*n*exp(-1i*k0*n*ytest);
        dy = 0;
    case 1i
        dx = 0;
        dy = -1i*k0*n*exp(1i*k0*n*xtest);
    case -1i    
        dx = 0;
        dy = 1i*k0*n*exp(-1i*k0*n*xtest);
end
    E = inc_Ez_field(source,test,n,k0,c,OMEGA,direction);    
        
 Hx = dx - E.*(shift_fact_x);       
 Hy = dy - E.*(shift_fact_y);       
       
end