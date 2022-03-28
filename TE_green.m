
function [G] = TE_green(source,test,n,k0,c,OMEGA, params)

a = source(1,:);
b = source(2,:);
x = test(1,:);
y = test(2,:);

% rot = exp(1i*k0*OMEGA/c*(a.*y-b.*x));
% Gst = 1i/4*besselh(0,1,k0*n*(sqrt((x-a).*(x-a)+(y-b).*(y-b))));

rr = sqrt((x-a)^2 + (y-b)^2);


gst = 1i/4*besselh(0,1,k0*n*sqrt(rr));
rot = exp(1i*k0*OMEGAc*(-b*x+a*y));

G_rot = gst.*rot;


const_1 = 1i/(params.omega*params.e0.*params.er_out);

etha = sqrt(params.mr_out/params.er_out);
const_2 = etha * params.OMEGA/params.c;

const_3 = 1i*params.omega*params.m0*params.mr_out*(params.OMEGA/params.c)^2;

%% calc first term

my_yy_1 = 1i/4*k0*n*besselh(1,1,k0*n*rr)*((2*(y-b)^2/(rr^3)-1/rr)+1i*k0*OMEGAc*(a+x)*(b-y)/rr)*rot;
my_yy_2 = 1i/4*(-k0^2*n^2*(b-y)^2/(rr^2)*besselh(0,1,k0*n*rr)+1i*k0*OMEGAc*x*1i*k0*OMEGAc*a*besselh(0,1,k0*n*rr))*rot;

A11 = my_yy_1 + my_yy_2; 

my_xx_1 = 1i/4*k0*n*besselh(1,1,k0*n*rr)*((2*(x-a)^2/(rr^3)-1/rr)+1i*k0*OMEGAc*(x-a)*(y+b)/rr)*rot;
my_xx_2 = 1i/4*(-k0^2*n^2*(a-x)^2/(rr^2)*besselh(0,1,k0*n*rr)+1i*k0*OMEGAc*y*1i*k0*OMEGAc*b*besselh(0,1,k0*n*rr))*rot;

A22 = my_xx_1 + my_xx_2;

my_xb_1 = -1i/4*k0*n*besselh(1,1,k0*n*rr)*(2*(b-y)*(a-x)/(rr^3)-(b-y)/rr*(1i*k0*OMEGAc*b)+1i*k0*OMEGAc*(a-x)/rr*x)*rot;
my_xb_2 = 1i/4*besselh(0,1,k0*n*rr)*((k0*n)^2*(b-y)*(a-x)/(rr^2)+(1i*k0*OMEGAc)^2*x*b-1i*k0*OMEGAc)*rot;

A21 = my_xb_1+ my_xb_2;

my_ya_1 = -1i/4*k0*n*besselh(1,1,k0*n*rr)*(2*(b-y)*(a-x)/(rr^3)+(a-x)/rr*(1i*k0*OMEGAc*a)-1i*k0*OMEGAc*(b-y)*y/rr)*rot;
my_ya_2 = 1i/4*besselh(0,1,k0*n*rr)*((k0*n)^2*(b-y)*(a-x)/(rr^2)+(1i*k0*OMEGAc)^2*y*a+1i*k0*OMEGAc)*rot;

A12 = my_ya_1 + my_ya_2;

FIRST_MAT =  [A11 , A12; A21, A22];

%%calc 2nd term





end

