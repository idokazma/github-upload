%% set up the sources location and test points locations. test points are on the surface of the scatterer;

function [E_SOL_st_fil_TM, E_sca] = filaments_TM(params,E_inc_z,H_inc_x,H_inc_y,Esym)


R = params.radius;
theta = linspace(0,2*pi,params.N_filaments+1);
theta = theta(1:end-1);
out_sources_loc = (params.R_out*R)*exp(1i*theta);
in_sources_loc =(params.R_in*R)*exp(1i*theta);

theta_test = linspace(0,2*pi,params.N_testpoints+1);
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
e0=params.e0;
mu0 = params.mu0;

mr_in = params.mr_in;
mr_out = params.mr_out;
er_in = params.er_in;
er_out = params.er_out;
lambda = params.lambda;


c=params.c;
f = params.f;
k0 = params.k0;
omega = params.omega;
OMEGA = params.OMEGA;

n_in = params.n_in;
n_out = params.n_out;

%% START THE ELECTROMAGNETIC CALCULATIONS

% set the Ez_incident and the Ht_incident

% now, for each Test point we are going to calculate the "impact" of a unit
% current located outside and inside the scatterer ON the surface (Test
% Points). as in Yehuda's article, we calculate first ZeI, then ZeII.
% ZeI (or ZeS) is the impact of the INSIDE currents, radiating with the
% OUTSIDE material parameters and Green function measured on the surface of
% the scatterer. ZeII is the OUTSIDE currents radiating towards the INSIDE
% with the INSIDE parameters.
syms xsym ysym zsym k0sym
Xsym = [xsym ysym zsym];
Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);



for i = 1:length(test_loc)
    E_inc(i,:) = double(subs(Esym,[xsym,ysym,k0sym] ,[test_loc(1,i),test_loc(2,i),k0]));
    H_inc(i,:) = double(subs(Hsym,[xsym,ysym,k0sym] ,[test_loc(1,i),test_loc(2,i),k0]));
end

H_inc_x = H_inc(:,1);
H_inc_y = H_inc(:,2);
H_inc_z = H_inc(:,3);

E_inc_x = E_inc(:,1);
E_inc_y = E_inc(:,2);
E_inc_z = E_inc(:,3);

for i=1:length(test_loc)
    %calculate the nX(E1-E2) scalar thing
    for j=1:length(in_sources_loc) % calculate the G function of the sources located *INSIDE* the scatterer at the test point;
        Mx(i,j) = 1i*omega*(mu0*mr_out)*scalar_green(in_sources_loc(:,j),test_loc(:,i),n_out,k0,c,OMEGA);
        
    end
    for k = 1:length(out_sources_loc)  % calculate the G function of the sources located *OUTSIDE* the scatterer at the test point, now with with the MINUS sign
        Mx(i,j+k) =  -1i*omega*(mu0*mr_in)*scalar_green(out_sources_loc(:,k),test_loc(:,i),n_in,k0,c,OMEGA);
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
    sol = cross(nhat(i,:),[H_inc_x(i),H_inc_y(i),0]);
    Vhz(i) = -1*sol(3);
end


M_total = [Mx;MHz];
figure; imagesc(real(M_total));

P_inverse_M = ((M_total')*M_total)\(M_total');
% P_inverse_M = (transpose(M_total)*M_total)\(transpose(M_total));

I_solution = P_inverse_M*(transpose([Vex,Vhz]));
% I_solution = (M_total)\(transpose([Vex,Vhz]));
% I_solution = linsolve(M_total,transpose([Vex,Vhz]));
%  I_solution =P_inverse_M*transpose([Vex,Vhz]);



% figure; quiver(test_loc(1,:),test_loc(2,:),direction(1,:),direction(2,:));
% figure; scatter(test_loc(1,:),test_loc(2,:),20,abs(E_inc_z),'filled');
% hold on;
% scatter(in_sources_loc(1,:),in_sources_loc(2,:),20,abs(I_solution(1:length(in_sources_loc))),'filled');
% scatter(out_sources_loc(1,:),out_sources_loc(2,:),20,abs(I_solution((1+length(in_sources_loc)):2*length(in_sources_loc))),'filled');
X = params.X;
Y = params.Y;
loc_line = [X(:),Y(:)]';
out_loc = (X.*X+Y.*Y)>R*R;
in_loc = (X.*X+Y.*Y)<=R*R;

out_loc_line = [X(out_loc==1),Y(out_loc==1)]';
in_loc_line = [X(in_loc==1),Y(in_loc==1)]';

E_sol = zeros(1,length(loc_line));

for tt=1:length(in_sources_loc)
    E_sol(out_loc(:)) = E_sol(out_loc(:))+1i*omega*(mu0*mr_out)*I_solution(tt)*scalar_green(in_sources_loc(:,tt),out_loc_line,n_out,k0,c,OMEGA);
end

for tt=1:length(out_sources_loc)
    E_sol(in_loc(:)) = E_sol(in_loc(:))+1i*omega*(mu0*mr_in)*I_solution(tt+length(in_sources_loc))*scalar_green(out_sources_loc(:,tt),in_loc_line,n_in,k0,c,OMEGA);
end

% tic
% inc_field = zeros(length(out_loc_line),3);
% for i=1:length(out_loc_line)
%     inc_field(i,:) =  double(subs(Esym,[xsym,ysym] ,[out_loc_line(1,i),out_loc_line(2,i)]));
%     
%     if mod(i,100)==0
%         est_t =toc;
%         (length(out_loc_line)-i)/100*est_t
%         tic
%     end
% end

Esym(3)

E_sol_total(out_loc(:)) = E_sol(out_loc(:)) +  exp(-1i*k0*out_loc_line(1,:));
E_sol_total(in_loc(:)) = E_sol(in_loc(:));


E_sol_sca(out_loc(:)) = E_sol(out_loc(:));
E_sol_sca(in_loc(:)) = E_sol(in_loc(:))-  exp(-1i*k0*in_loc_line(1,:));
%
% figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol),'filled'); axis equal
% figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol),'filled');axis equal
%
%
% figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol_maybe),'filled'); axis equal
% figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol_maybe),'filled');axis equal

% figure; imagesc(abs(reshape(E_sol_maybe,length(X),length(X)))); axis equal; title('filaments -- ABS')
% 
% figure; imagesc(imag(reshape(E_sol_maybe,length(X),length(X)))); axis equal ; title('filaments --IMAG')
% figure; imagesc(real(reshape(E_sol_maybe,length(X),length(X)))); axis equal ;  title('filaments --REAL')
E_SOL_st_fil_TM = reshape(E_sol_total,length(X),length(X));
E_sca = reshape(E_sol_sca,length(X),length(X));
end
