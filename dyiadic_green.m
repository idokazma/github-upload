function [Gx,Gy] = dyiadic_green(source,test,params, inside)

if (inside)
n = params.n_in;
mu = (params.mu0*params.mr_in);
else
n = params.n_out;
mu = (params.mu0*params.mr_out); 
end

c = params.c;
OMEGA = params.OMEGA;
k0 = params.k0;
Iz = params.Iz;
omega = params.omega;



xsource = source(1,:);
ysource = source(2,:);
xtest = test(1,:);
ytest = test(2,:);

rot = exp(1i*k0*OMEGA/c*(xsource.*ytest-ysource.*xtest));

% rot_derx = +1i*k0*OMEGA/c.*xsource.*rot;
% rot_dery = +1i*k0*OMEGA/c.*ysource.*rot;
% 
rr =sqrt((xtest-xsource).*(xtest-xsource)+(ytest-ysource).*(ytest-ysource));
% 
% Gst = 1i/4*besselh(0,1,k0*n*(rr));
% 
% Gx = 1i/4*(-1*k0*n*(ytest-ysource)./(rr)).*besselh(1,1,k0*n*rr);
% Gy = 1i/4*( 1*k0*n*(xtest-xsource)./(rr)).*besselh(1,1,k0*n*rr);
% 
% shift_fact_x = -1i*k0*OMEGA/c.*xtest;
% shift_fact_y = -1i*k0*OMEGA/c.*ytest;
% 
% Gx_rot = Gx.*rot + Gst.*rot_derx + shift_fact_x.*Gst.*rot;
% Gy_rot = Gy.*rot + Gst.*rot_dery + shift_fact_y.*Gst.*rot;
% 
% fact = 1/(4*1i);
% Gx_1 = (k0*n*(ytest-ysource)./(rr)).*besselh(1,1,k0*n*rr);
% Gx_2 = -k0^2*OMEGA/omega*(xsource-xtest).*besselh(0,1,k0*n*rr);
% 
% Gy_1 = (k0*n*(xtest-xsource)./(rr)).*besselh(1,1,k0*n*rr);
% Gy_2 = k0^2*OMEGA/omega*(ysource-ytest).*besselh(0,1,k0*n*rr);


fact_1 = 1i/4*k0*n;
fact_2 = OMEGA*omega/(4*c*c);

Gx_1 = ((ytest-ysource)./(rr)).*(-1).*besselh(1,1,k0*n*rr);
Gx_2 = (xtest-xsource).*besselh(0,1,k0*n*rr);

Gy_1 = (-(xtest-xsource)./(rr)).*(-1).*besselh(1,1,k0*n*rr);
Gy_2 = (ytest-ysource).*besselh(0,1,k0*n*rr);

Gx = fact_1 * Gx_1 + fact_2 * Gx_2;
Gy = fact_1 * Gy_1 + fact_2 * Gy_2;

Gx = Gx .* rot;
Gy = Gy .* rot;


end