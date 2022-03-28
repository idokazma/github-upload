clear all
%% set up the sources location and test points locations. test points are on the surface of the scatterer;
R = 0.001;
theta = linspace(0,2*pi,81);
theta = theta(1:end-1);
out_sources_loc = (1.2*R)*exp(1i*theta);
in_sources_loc =(0.8*R)*exp(1i*theta);

theta_test = linspace(0,2*pi,501);
theta_test = theta_test(1:end-1);
test_loc = (1*R)*exp(1i*theta_test);

%Transform to X,Y
out_sources_loc = [real(out_sources_loc);imag(out_sources_loc)];
in_sources_loc = [real(in_sources_loc);imag(in_sources_loc)];
test_loc = [real(test_loc);imag(test_loc)];

%calc the direction vector of the surface border (nhat)
for i = 1:length(test_loc)
    direction(:,i) = test_loc(:,i)/norm(test_loc(:,i));
    nhat(i,:) = [direction(1,i),direction(2,i),0];
end


%% set some EM parameters of the area
e0=8.85418781762039*1e-12;
mu0 = 4*pi*1e-7;

mr_in = 10;
mr_out = 1;
er_in = 1;
er_out = 1;

c=1/sqrt(e0*mu0);
lambda = 1;
f = c/lambda;
k0 = 2*pi/lambda;
omega = f*2*pi;

n_in = sqrt(mr_in*er_in);
n_out = sqrt(mr_out*er_out);

%% START THE ELECTROMAGNETIC CALCULATIONS

% set the Ez_incident and the Ht_incident
E_inc_z = exp(-1i*k0*test_loc(1,:));
%  H_inc_y = -sqrt(e0/mu0)*exp(-1i*k0*test_loc(1,:));
 H_inc_y =1i/(omega*(mu0*mr_out))*(-1i*k0*exp(-1i*k0*test_loc(1,:)));

H_inc_x = 0*(H_inc_y);

% now, for each Test point we are going to calculate the "impact" of a unit
% current located outside and inside the scatterer ON the surface (Test
% Points). as in Yehuda's article, we calculate first ZeI, then ZeII.
% ZeI (or ZeS) is the impact of the INSIDE currents, radiating with the
% OUTSIDE material parameters and Green function measured on the surface of
% the scatterer. ZeII is the OUTSIDE currents radiating towards the INSIDE
% with the INSIDE parameters.

for i=1:length(test_loc)
     %calculate the nX(E1-E2) scalar thing
        for j=1:length(in_sources_loc) % calculate the G function of the sources located *INSIDE* the scatterer at the test point;
            Mx(i,j) = 1i*omega*(mu0*mr_out)*scalar_green(in_sources_loc(:,j),test_loc(:,i),n_out,k0);

        end
        for k = 1:length(out_sources_loc)  % calculate the G function of the sources located *OUTSIDE* the scatterer at the test point, now with with the MINUS sign
            Mx(i,j+k) =  -1i*omega*(mu0*mr_in)*scalar_green(out_sources_loc(:,k),test_loc(:,i),n_in,k0);
        end
        Vex(i) = -E_inc_z(i);
    
end


for i=1:length(test_loc) %% calculate nXdeltaH(Z direction)
    for j=1:length(in_sources_loc)
        [Gx,Gy] = dyiadic_green(in_sources_loc(:,j),test_loc(:,i),n_out,k0);
        sol = cross(nhat(i,:),[Gx,Gy,0]);
        MHz(i,j) = sol(3);
%         MHz(i,j) = j;
    end
    
    for k=1:length(out_sources_loc)
        [Gx,Gy] = dyiadic_green(out_sources_loc(:,k),test_loc(:,i),n_in,k0);
        sol = cross(nhat(i,:),[Gx,Gy,0]);
        MHz(i,j+k) = -1*sol(3);
%         MHz(i,j+k) = -k;
    end
    sol = cross(nhat(i,:),[0,H_inc_y(i),0]);
    Vhz(i) = -1*sol(3);
end


M_total = [Mx;MHz];
P_inverse_M = ((M_total')*M_total)\(M_total');
% P_inverse_M = (transpose(M_total)*M_total)\(transpose(M_total));

I_solution = P_inverse_M*(transpose([Vex,Vhz]));
% I_solution = (M_total)\(transpose([Vex,Vhz]));
% I_solution = linsolve(M_total,transpose([Vex,Vhz]));
%  I_solution =P_inverse_M*transpose([Vex,Vhz]);



figure; quiver(test_loc(1,:),test_loc(2,:),direction(1,:),direction(2,:));
figure; scatter(test_loc(1,:),test_loc(2,:),20,abs(E_inc_z),'filled');
hold on;
scatter(in_sources_loc(1,:),in_sources_loc(2,:),20,abs(I_solution(1:length(in_sources_loc))),'filled');
scatter(out_sources_loc(1,:),out_sources_loc(2,:),20,abs(I_solution((1+length(in_sources_loc)):2*length(in_sources_loc))),'filled');

len_n = 480;
wid_n = 480;

len = 0.01;
wid = 0.01;

y = linspace(0,len,len_n)- len/2;
x = linspace(0,wid,wid_n)- wid/2;
[X,Y] = meshgrid(x,y);
loc_line = [X(:),Y(:)]';
out_loc = (X.*X+Y.*Y)>R*R;
in_loc = (X.*X+Y.*Y)<=R*R;

out_loc_line = [X(out_loc==1),Y(out_loc==1)]';
in_loc_line = [X(in_loc==1),Y(in_loc==1)]';

E_sol = zeros(1,length(loc_line));

for tt=1:length(in_sources_loc)
    E_sol(out_loc(:)) = E_sol(out_loc(:))+1i*omega*(mu0*mr_out)*I_solution(tt)*scalar_green(in_sources_loc(:,tt),out_loc_line,n_out,k0);
end

for tt=1:length(out_sources_loc)
    E_sol(in_loc(:)) = E_sol(in_loc(:))+1i*omega*(mu0*mr_in)*I_solution(tt+length(in_sources_loc))*scalar_green(out_sources_loc(:,tt),in_loc_line,n_in,k0);
end

E_sol_maybe(out_loc(:)) = E_sol(out_loc(:)) + exp(-1i*k0*out_loc_line(1,:));
E_sol_maybe(in_loc(:)) = E_sol(in_loc(:));

% 
% figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol),'filled'); axis equal
% figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol),'filled');axis equal
% 
% 
% figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol_maybe),'filled'); axis equal
% figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol_maybe),'filled');axis equal

figure; imagesc(abs(reshape(E_sol_maybe,length(X),length(X)))); axis equal; title('filaments -- ABS')

figure; imagesc(imag(reshape(E_sol_maybe,length(X),length(X)))); axis equal ; title('filaments --IMAG')
figure; imagesc(real(reshape(E_sol_maybe,length(X),length(X)))); axis equal ;  title('filaments --REAL')

