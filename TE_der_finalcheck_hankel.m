syms x y z a b c k0 n OMEGAc
%f = besselj(1,kn*sqrt((x-a)^2 + (y-b)^2  ));
% f = (b-y)/sqrt((x-a)^2 + (y-b)^2  );
% f = besselj(1,kn*sqrt((x-a)^2 + (y-b)^2  ));
% pretty((diff(f,y)))
gst = 1i/4*besselh(0,1,k0*n*sqrt((x-a)^2 + (y-b)^2  ));
rot = exp(1i*k0*OMEGAc*(-b*x+a*y));

V_test = [x,y,0];
V_source = [a,b,0];

g_rot = gst*rot;

pretty(gst)
pretty(rot)
pretty(g_rot)

disp ('~~~~~~~~~~~~Start dx~~~~~~~~~~~~~~');
resx = diff(g_rot,x);
compx = -1i/4*(k0*n*(x-a)/(sqrt((x-a)^2 + (y-b)^2  ))*besselh(1,1,k0*n*sqrt((x-a)^2 + (y-b)^2))+1i*k0*OMEGAc*b*besselh(0,1,k0*n*sqrt((x-a)^2 + (y-b)^2  )))*rot;
disp ('dx - My result');
pretty(compx)
disp ('dx - MATLAB result');
pretty(resx)
disp ('Me vs MATLAB (0 means IDENTICAL)');
pretty(simplify(compx-resx))
disp ('~~~~~~~~~~~~~~END dx~~~~~~~~~~~~~');

disp ('~~~~~~~~~~~~Start dy~~~~~~~~~~~~~~');
resy = diff(g_rot,y);
compy = -1i/4*(k0*n*(y-b)/(sqrt((x-a)^2 + (y-b)^2  ))*besselh(1,1,k0*n*sqrt((x-a)^2 + (y-b)^2  ))-1i*k0*OMEGAc*a*besselh(0,1,k0*n*sqrt((x-a)^2 + (y-b)^2  )))*rot;
% disp ('dy - My result');
% pretty(compy)
% disp ('dy - MATLAB result');
% pretty(resy)
disp ('Me vs MATLAB (0 means IDENTICAL)');
pretty(simplify(compy-resy))
disp ('~~~~~~~~~~~~~~END dy~~~~~~~~~~~~~');

disp ('~~~~~~~~~~~~Start da~~~~~~~~~~~~~~');
resa = diff(g_rot,a);
compa = -1i/4*(k0*n*(a-x)/(sqrt((x-a)^2 + (y-b)^2  ))*besselh(1,1,k0*n*sqrt((x-a)^2 + (y-b)^2  ))-1i*k0*OMEGAc*y*besselh(0,1,k0*n*sqrt((x-a)^2 + (y-b)^2  )))*rot;
% disp ('dy - My result');
% pretty(compa)
% disp ('dy - MATLAB result');
% pretty(resa)
disp ('Me vs MATLAB (0 means IDENTICAL)');
pretty(simplify(compa-resa))
disp ('~~~~~~~~~~~~~~END da~~~~~~~~~~~~~');


disp ('~~~~~~~~~~~~Start db~~~~~~~~~~~~~~');
resb = diff(g_rot,b);
compb = -1i/4*(k0*n*(b-y)/(sqrt((x-a)^2 + (y-b)^2  ))*besselh(1,1,k0*n*sqrt((x-a)^2 + (y-b)^2  ))+1i*k0*OMEGAc*x*besselh(0,1,k0*n*sqrt((x-a)^2 + (y-b)^2  )))*rot;
% disp ('dy - My result');
% pretty(compb)
% disp ('dy - MATLAB result');
% pretty(resb)
disp ('Me vs MATLAB (0 means IDENTICAL)');
pretty(simplify(compb-resb))
disp ('~~~~~~~~~~~~~~END db~~~~~~~~~~~~~');

% simplify(compx-resx)
% simplify(compy-resy)
% simplify(compa-resa)
% simplify(compb-resb)


rr = sqrt((x-a)^2 + (y-b)^2);


% my_xx_1 = 1i/4*k0*n*besselj(1,k0*n*rr)*((2*(x-a)^2/(rr^3)-1/rr)+1i*k0*OMEGAc*(x-a)*(y+b)/rr)*rot;
% my_xx_2 = 1i/4*(-k0^2*n^2*(a-x)^2/(rr^2)*besselj(0,k0*n*rr)+1i*k0*OMEGAc*y*1i*k0*OMEGAc*b*besselj(0,k0*n*rr))*rot;
% 
% pretty(simplify((my_xx_1+my_xx_2)-simplify(-diff(diff(g_rot,a),x))))
% 
% disp ('~~~~~~~~~~~~dyy~~~~~~~~~~~~~~');
% my_yy_1 = 1i/4*k0*n*besselj(1,k0*n*rr)*((2*(y-b)^2/(rr^3)-1/rr)+1i*k0*OMEGAc*(a+x)*(b-y)/rr)*rot;
% my_yy_2 = 1i/4*(-k0^2*n^2*(b-y)^2/(rr^2)*besselj(0,k0*n*rr)+1i*k0*OMEGAc*x*1i*k0*OMEGAc*a*besselj(0,k0*n*rr))*rot;
% 
% pretty(simplify((my_yy_1+my_yy_2)-simplify(-diff(diff(g_rot,b),y))))
% 
% disp ('~~~~~~~~~~~~dxy^~~~~~~~~~~~~~~');
% my_xb_1 = -1i/4*k0*n*besselj(1,k0*n*rr)*(2*(b-y)*(a-x)/(rr^3)-(b-y)/rr*(1i*k0*OMEGAc*b))*rot-1i*(-1i/4)*k0*OMEGAc*k0*n*(x-a)*x/rr*besselj(1,k0*n*rr)*rot;
% my_xb_2 = 1i/4*k0*n*(k0*n*(b-y)*(a-x)/(rr^2)*besselj(0,k0*n*rr))*rot+1i/4*(1i*k0*OMEGAc*x*1i*k0*OMEGAc*b*besselj(0,k0*n*rr)*rot-1i*k0*OMEGAc*besselj(0,k0*n*rr)*rot);
% 
% pretty((simplify((my_xb_1+my_xb_2)-simplify(diff(diff(g_rot,b),x)))))
% 
% disp ('~~~~~~~~~~~~dyx^~~~~~~~~~~~~~~');
% 
% my_ya_1 = -1i/4*k0*n*besselj(1,k0*n*rr)*(2*(b-y)*(a-x)/(rr^3)+(a-x)/rr*(1i*k0*OMEGAc*a))*rot+1i*(-1i/4)*k0*OMEGAc*k0*n*(y-b)*y/rr*besselj(1,k0*n*rr)*rot;
% my_ya_2 = 1i/4*k0*n*(k0*n*(b-y)*(a-x)/(rr^2)*besselj(0,k0*n*rr))*rot+1i/4*(1i*k0*OMEGAc*y*1i*k0*OMEGAc*a*besselj(0,k0*n*rr)*rot+1i*k0*OMEGAc*besselj(0,k0*n*rr)*rot);
% 
% pretty((simplify((my_ya_1+my_ya_2)-simplify(diff(diff(g_rot,a),y)))))
% 


disp ('~~~~~~~~~~~~dxx~~~~~~~~~~~~~~');

my_xx_1 = 1i/4*k0*n*besselh(1,1,k0*n*rr)*((2*(x-a)^2/(rr^3)-1/rr)+1i*k0*OMEGAc*(x-a)*(y+b)/rr)*rot;
my_xx_2 = 1i/4*(-k0^2*n^2*(a-x)^2/(rr^2)*besselh(0,1,k0*n*rr)+1i*k0*OMEGAc*y*1i*k0*OMEGAc*b*besselh(0,1,k0*n*rr))*rot;

pretty(simplify((my_xx_1+my_xx_2)-simplify(-diff(diff(g_rot,a),x))))


disp ('~~~~~~~~~~~~dyy~~~~~~~~~~~~~~');
my_yy_1 = 1i/4*k0*n*besselh(1,1,k0*n*rr)*((2*(y-b)^2/(rr^3)-1/rr)+1i*k0*OMEGAc*(a+x)*(b-y)/rr)*rot;
my_yy_2 = 1i/4*(-k0^2*n^2*(b-y)^2/(rr^2)*besselh(0,1,k0*n*rr)+1i*k0*OMEGAc*x*1i*k0*OMEGAc*a*besselh(0,1,k0*n*rr))*rot;

pretty(simplify((my_yy_1+my_yy_2)-simplify(-diff(diff(g_rot,b),y))))
disp ('~~~~~~~~~~~~dxy^~~~~~~~~~~~~~~');
my_xb_1 = -1i/4*k0*n*besselh(1,1,k0*n*rr)*(2*(b-y)*(a-x)/(rr^3)-(b-y)/rr*(1i*k0*OMEGAc*b)+1i*k0*OMEGAc*(a-x)/rr*x)*rot;
my_xb_2 = 1i/4*besselh(0,1,k0*n*rr)*((k0*n)^2*(b-y)*(a-x)/(rr^2)+(1i*k0*OMEGAc)^2*x*b-1i*k0*OMEGAc)*rot;

pretty((simplify((my_xb_1+my_xb_2)-simplify(diff(diff(g_rot,b),x)))))
disp ('~~~~~~~~~~~~dyx^~~~~~~~~~~~~~~');
my_ya_1 = -1i/4*k0*n*besselh(1,1,k0*n*rr)*(2*(b-y)*(a-x)/(rr^3)+(a-x)/rr*(1i*k0*OMEGAc*a)-1i*k0*OMEGAc*(b-y)*y/rr)*rot;
my_ya_2 = 1i/4*besselh(0,1,k0*n*rr)*((k0*n)^2*(b-y)*(a-x)/(rr^2)+(1i*k0*OMEGAc)^2*y*a+1i*k0*OMEGAc)*rot;

pretty((simplify((my_ya_1+my_ya_2)-simplify(diff(diff(g_rot,a),y)))))
