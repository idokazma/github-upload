%% Mie Series Calculator
function [E_in, E_SOL_st_mie,current_mie] = Mie_Series(params)


len_n = params.len_n;
wid_n = params.wid_n;

len = params.len;
wid = params.wid;

y = linspace(0,len,len_n)- len/2;
x = linspace(0,wid,wid_n)- wid/2;


[X,Y] = meshgrid(x,y);

total_cells = len_n*wid_n;
epsilon_bg=1;
lambda = params.lambda;
k0 = 2*pi/lambda;

beta = 1;

V_i = abs(x(2)-x(1))*abs(y(2)-y(1));
R_i_eff = sqrt(V_i/pi); %% R_i_eff = sqrt(V_i/pi);
gamma=R_i_eff/k0*besselh(1,1,R_i_eff*k0)+1i*2/(pi*k0^2);

%% source gen
E_bg_dy = zeros(size(X));
E_bg_st = zeros(size(X));

% source_location = [5*ones(1,15);linspace(3,-3,15)];%;72.2,76.3];
source_location = [5000000;0];%;72.2,76.3];

c = 3*10^8;

omega = c/lambda;
tic
% Omega_B=omega*0.005;
% for(ll=1:size(source_location,2))
%  
%     
% %     E_bg_dy = E_bg_dy + 1i/4*besselh(0,1,k0.*sqrt((X-source_location(1,ll)).^2+(Y-source_location(2,ll)).^2)).*exp(1i*Omega_B/omega*k0^2.*(X*source_location(2,ll)-Y*source_location(1,ll)));
%     E_bg_st = E_bg_st + 1i/4*besselh(0,1,k0.*sqrt((X-source_location(1,ll)).^2+(Y-source_location(2,ll)).^2));
% 
% end
% toc
% E_bg_st = E_bg_st/abs(max(max(E_bg_st)));

E_bg_st = exp(-1i*k0*X);

% E_bg_st = E_bg_st*exp(-1i*phase(E_bg_st(125,125)));
% figure; imagesc(real(E_bg_dy)); set(gca,'YDir','normal'); axis equal;
% figure; imagesc(real(E_bg_st)); set(gca,'YDir','normal'); axis equal;
% 
% figure; imagesc(abs(E_bg_st).^2); set(gca,'YDir','normal'); axis equal;

scattereres = epsilon_bg*ones(wid_n,len_n);
% 
% center = [0,0];
% radius = 5.2;
% loc =find(and(((X-center(1)).^2 + (Y-center(2)).^2)<radius^2,((X-center(1)).^2 + (Y-center(2)).^2)>(radius-0.4)^2));

 center = [0,0];
 radius = params.radius;
 loc =find(sqrt((X-center(1)).^2 + (Y-center(2)).^2)<radius);

 mr_in = params.mr_in;
mr_out = params.mr_out;
er_in = params.er_in;
er_out = params.er_out;
 

delta_eps = params.delta_eps;
% scattereres (3,3) = 2;
% scattereres (3,500) = 2;

[theta,rho] = cart2pol(X,Y);

E_final=zeros(size(X));
% E_0 = abs(max(max(E_bg_st)));
E_0 = (E_bg_st(x==0,y==0));

n_1 = sqrt(er_out*mr_out);
n_2 = sqrt(er_in*mr_in);

mu=1;
etha_1 = sqrt(mr_out/er_out);
etha_2 = sqrt(mr_in/er_in);

m=0:50;
tic
%% b_tm calc
b_m_tm_up_A = etha_2*(1/2*(besselj(m-1,k0*n_1*radius)-besselj(m+1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_tm_up_B = etha_1*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselj(m,k0*n_1*radius);

b_m_tm_down_C = etha_2*(1/2*(besselh(m-1,1,k0*n_1*radius)-besselh(m+1,1,k0*n_1*radius))).*besselj(m,k0*n_2*radius);
b_m_tm_down_D = etha_1*(1/2*(besselj(m-1,k0*n_2*radius)-besselj(m+1,k0*n_2*radius))).*besselh(m,1,k0*n_1*radius);

b_TM = ((-1)*(-1i).^m).*(b_m_tm_up_A-b_m_tm_up_B)./(b_m_tm_down_C-b_m_tm_down_D);

%% c_tm calc

c_TM = b_TM.*besselh(m,1,k0*n_1*radius)./besselj(m,k0*n_2*radius)+(-1i).^m.*besselj(m,k0*n_1*radius)./besselj(m,k0*n_2*radius);


sigma =ones(size(m))*2;
sigma(1)=1;
E_mid=E_final;
for(i=1:length(m))
E_mid(rho>radius) =  E_mid(rho>radius)+sigma(i)*b_TM((i)).*besselh(m(i),1,k0*n_1.*rho(rho>radius)).*cos(m(i)*theta(rho>radius));
 E_mid(rho<=radius) =  E_mid(rho<=radius)+sigma(i)*c_TM(i).*besselj(m(i),k0*n_2.*rho(rho<=radius)).*cos(m(i)*theta(rho<=radius));

end
E_mid(rho>radius) = E_0*E_mid(rho>radius)+(E_bg_st(rho>radius));
E_mid(rho<=radius)= E_0*E_mid(rho<=radius);
toc
% figure; imagesc(real(E_mid)); axis equal
% figure; imagesc(x,y,abs(E_mid).^2,[0 10]); set(gca,'YDir','normal'); axis equal; title ('ENG Mie Series');
% figure; imagesc(db((E_mid).^2)); set(gca,'YDir','normal'); axis equal


% figure; imagesc(db(E_mid)); axis equal; (find(max(max(db(E_mid)))==db(E_mid)))
% figure; imagesc(db(E_bg_st)); axis equal

% figure; subplot(1,2,1); imagesc(x,y,abs(E_mid).^2,[0 4]); set(gca,'YDir','normal'); axis equal; title ('ENG Mie Series');  axis image; xlabel('x(\lambda)'); ylabel('y(\lambda)'); colorbar
% subplot(1,2,2); imagesc(x,y,abs(SOL_st).^2, [0 4]); set(gca,'YDir','normal'); axis equal; title ('ENG MOM');  axis image; xlabel('x(\lambda)'); ylabel('y(\lambda)'); colorbar
% 
% loc_x = find(x==0);
% loc_y = find(y==0);
% figure; plot(Y(:,loc_x),real(E_mid(:,loc_x)),'.-'); hold on; plot(Y(:,loc_x),real(SOL_st(:,loc_x)),'.-');
% figure; plot(Y(:,loc_x),imag(E_mid(:,loc_x)),'.-'); hold on; plot(Y(:,loc_x),imag(SOL_st(:,loc_x)),'.-');
% 
% figure; plot(X(loc_y,:),real(E_mid(loc_y,:))); hold on; plot(X(loc_y,:),real(SOL_st(loc_y,:)));
% figure; plot(X(loc_y,:),imag(E_mid(loc_y,:))); hold on; plot(X(loc_y,:),imag(SOL_st(loc_y,:)));

figure; imagesc(real(E_mid)); axis equal; title('Mie - REAL');
figure; imagesc(imag(E_mid)); axis equal; title('Mie - IMAG');
figure; imagesc(abs(E_mid)); axis equal; title('Mie - ABS');


E_in = E_bg_st;
E_SOL_st_mie = E_mid;
current_mie = sum(E_mid(loc))*V_i;