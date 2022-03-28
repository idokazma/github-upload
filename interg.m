R = 1/100*(-1:1/10:1)*params.lambda;
theta = linspace(0,pi,25);
theta = [theta];
Rc = 10000.001*params.lambda;
res = zeros (length(R), length(theta),length(R), length(theta));
bes = res;

f = waitbar(1,'wish for good');
for a = 1:length(R)
    
    for b = 1:length(theta)
        
        
        for c = 1:length(R)
            
            for d = 1:length(theta)
                
                 A = [Rc + R(a)*cos(theta(b)) R(a)*sin(theta(b)) 0];
                  B = [1.000001*Rc + R(c)*cos(theta(d)) R(c)*sin(theta(d)) 0];
%                 B = [Rc + R(c)*cos(theta(d)) R(c)*sin(theta(d)) 0];
%                 C = cross(A,B);
                res(a,b,c,d) = (Rc + R(a)*cos(theta(b)))*(R(c)*sin(theta(d))) -(R(a)*sin(theta(b)))*(Rc + R(c)*cos(theta(d))) ;
                 bes(a,b,c,d) = 1i/4*besselh(0,1,params.k0*params.n_in*norm(A - B));
            end
            
        end
        
        
    end
    
    waitbar(1-(a/length(R))^4,f,'Finishing');
    
end
close (f);
Rc
inte = res.*res.*bes;
% disp('geo')
% mean(mean(mean(mean(res.*res))))
% disp('bes')
% mean(mean(mean(mean(bes))))
% disp('inte')
% mean(mean(mean(mean(inte))))

fprintf('geometric %d ; hankel %d ; integral %d ; \n', mean(mean(mean(mean(res.*res)))), mean(mean(mean(mean(bes)))),mean(mean(mean(mean(inte)))))

