
function [E] = inc_Ez_field(source,test,n,k0,c,OMEGA,direction)

% xsource = source(1,:);
% ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

switch direction
    
    case 1
        E = exp(1i*k0*n*ytest);
    case -1    
        E = exp(-1i*k0*n*ytest);
    case 1i
        E = exp(1i*k0*n*xtest);
    case -1i    
        E = exp(-1i*k0*n*xtest);
end   
 E = E;

end