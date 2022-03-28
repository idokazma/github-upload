%% set up the sources location and test points locations. test points are on the surface of the scatterer;

function [H_SOL_st_fil_TE,H_sol_sca, Ex_sol_final, Ey_sol_final, Ex_sol_final_sca, Ey_sol_final_sca ,alpha_mat ,output_mat] = filaments_TE_multiple(params,Hsym)

 H_SOL_st_fil_TE = [];
 H_sol_sca = [];
 Ex_sol_final = [];
 Ey_sol_final = [];
 Ex_sol_final_sca = [];
 Ey_sol_final_sca = [];
 alpha_mat  = [];


filaments_total = [];
test_point_total = [];
for k = 1 : length(params.sca_x)
    R = params.radius;
    theta = linspace(0,2*pi,params.N_filaments+1);
    theta = theta(1:end-1);
    out_sources_loc = (params.R_out*R)*exp(1i*theta);
    in_sources_loc =(params.R_in*R)*exp(1i*theta);
    
    theta_test = linspace(0,2*pi,params.N_testpoints+1);
    theta_test = theta_test(1:end-1);
    test_loc = (1*R)*exp(1i*theta_test);
    
    %Transform to X,Y
    
    out_sources_loc = [real(out_sources_loc)+params.sca_x(k);imag(out_sources_loc)+params.sca_y(k)];
    in_sources_loc = [real(in_sources_loc)+params.sca_x(k);imag(in_sources_loc)+params.sca_y(k)];
    
    filaments_out = [out_sources_loc;ones(1,length(out_sources_loc))*k;zeros(1,length(out_sources_loc))];
    filaments_in = [in_sources_loc;ones(1,length(in_sources_loc))*k;ones(1,length(in_sources_loc))];
    filaments_total = [filaments_total , filaments_in , filaments_out];
    
    test_loc = [real(test_loc);imag(test_loc)];
    
    %calc the direction vector of the surface border (nhat)
    for i = 1:length(test_loc)
        direction(:,i) = test_loc(:,i)/norm(test_loc(:,i));
        nhat(i,:) = [direction(1,i),direction(2,i),0];
    end
    test_loc = test_loc+[params.sca_x(k);params.sca_y(k)];
    
    test_point = [test_loc ; ones(1,length(test_loc))*k; nhat'];
    
    test_point_total = [test_point_total, test_point];
    
    
end
filaments_total = transpose(filaments_total);
test_point_total = transpose(test_point_total);

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

% Esym = 1i/(omega*(e0*er_out))*curl(Hsym,Xsym);

syms xsym ysym zsym k0sym
Xsym = [xsym ysym zsym];
if (false)
    
    for i = 1:length(test_point_total)
        E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),k0]));
        H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
        
        %     E_inc_green = scalar_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0,params.c,params.OMEGA);
        %     [H_inc_green_x , H_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0, c,params.OMEGA);
        %
    end
    disp ('End Eval Symbolic');
    toc
end


if (params.is_plane_wave == 1)
    
    
    for i = 1:length(test_point_total)
        E_inc(i,:) = double(subs(params.Esym,[xsym,ysym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),k0]));
        H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
    end
    
    H_inc_x = H_inc(:,1);
    H_inc_y = H_inc(:,2);
    H_inc_z = H_inc(:,3);
    
    E_inc_x = E_inc(:,1);
    E_inc_y = E_inc(:,2);
    E_inc_z = E_inc(:,3);
    
else
    
    H_inc_green = scalar_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0,params.c,params.OMEGA);
    [E_inc_green_x , E_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0, c,params.OMEGA);
    E_inc_green_x = -1 * E_inc_green_x*1./(1i*params.omega*(params.e0*params.er_out));
    E_inc_green_y = -1 * E_inc_green_y*1./(1i*params.omega*(params.e0*params.er_out));
    
    H_inc_z = H_inc_green;
    H_inc_x = 0*H_inc_z;
    H_inc_y = 0*H_inc_z;
    
    E_inc_x = E_inc_green_x;
    E_inc_y = E_inc_green_y;
    E_inc_z = 0*E_inc_green_x;
    
    
end

% %% just for test
%     H_inc_z = H_inc_green*1 + 0*mean(H_inc_green);
%     H_inc_x = 0*H_inc_green;
%     H_inc_y = 0*H_inc_green;
%     
%     E_inc_x = E_inc_green_x*1 + 0*mean(E_inc_green_x);
%     E_inc_y = E_inc_green_y*1 + 0*mean(E_inc_green_y);
%     E_inc_z = 0*E_inc_green_x;

% %      H_inc_green = scalar_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0,params.c,params.OMEGA);
% %     [E_inc_green_x , E_inc_green_y] = dyiadic_green([params.source_loc_x;params.source_loc_y],test_point_total',params.n_out,params.k0, c,params.OMEGA);
% %     E_inc_green_x = -1 * E_inc_green_x*1./(1i*params.omega*(params.e0*params.er_out));
% %     E_inc_green_y = -1 * E_inc_green_y*1./(1i*params.omega*(params.e0*params.er_out));
% %     
%     H_inc_z = 1+0*H_inc_green;
%     H_inc_x = 0*H_inc_z;
%     H_inc_y = 0*H_inc_z;
%     
%     E_inc_x = 0*E_inc_green_x;
%     E_inc_y = 1./(1i*params.omega*(params.e0*params.er_out))*(-1i*params.k0) + 0*E_inc_green_y;
%     E_inc_z = 0*E_inc_green_x;


for i=1:length(test_point_total)
    nhat = test_point_total(i,4:6);
    test_point_xy = [test_point_total(i,1:2)];
    
    for j=1:length(filaments_total)
        if (test_point_total(i,3) == filaments_total(j,3) && filaments_total(j,4) == 0) %% we are talking on the same scatterer and test point, and the filaments is outside the sca
            MHz(i,j) = -1i*omega*(e0*er_in)*scalar_green(filaments_total(j,1:2)',test_point_xy',n_in,k0,c,OMEGA);
            [Gx,Gy] =  dyiadic_green(filaments_total(j,1:2)',test_point_xy',n_in,k0,c,params.OMEGA);
            Gx = -Gx;
            Gy = -Gy;
            sol = cross(nhat,[Gx,Gy,0]);
            MEt(i,j) = -sol(3);%MINUS SIGN
            
        elseif (test_point_total(i,3) == filaments_total(j,3) && filaments_total(j,4) == 1)  %% we are talking on the same scatterer and test point, and the filaments is inside the sca
            MHz(i,j) = 1i * omega * (e0*er_out) * scalar_green(filaments_total(j,1:2)',test_point_xy',n_out,k0,c,OMEGA);
            [Gx,Gy] =  dyiadic_green(filaments_total(j,1:2)',test_point_xy',n_out,k0, c,params.OMEGA);
            Gx = -Gx;
            Gy = -Gy;
            sol = cross(nhat,[Gx,Gy,0]);
            
            MEt(i,j) = sol(3);
            
        elseif (test_point_total(i,3) ~= filaments_total(j,3) && filaments_total(j,4) == 1)  %% we are talking on other scatterer and test point, and the filaments is inside the sca
            MHz(i,j) = 1i * omega * (e0*er_out) * scalar_green(filaments_total(j,1:2)',test_point_xy',n_out,k0,c,OMEGA);
            [Gx,Gy] = dyiadic_green(filaments_total(j,1:2)',test_point_xy',n_out,k0,c,params.OMEGA);
            Gx = -Gx;
            Gy = -Gy;
            sol = cross(nhat,[Gx,Gy,0]);
            MEt(i,j) = sol(3);
            
        else
            MHz(i,j) = 0;
            MEt(i,j) = 0;
            
            
        end
        
    end
    Vhz(i) = -H_inc_z(i);
    sol = cross(nhat,[E_inc_x(i),E_inc_y(i),0]);
    Vet(i) = -1*sol(3);
end

% % % % % % % %
% % % % % % % %
% % % % % % % % for i=1:length(test_loc)
% % % % % % % %     %calculate the nX(E1-E2) scalar thing
% % % % % % % %     for j=1:length(in_sources_loc) % calculate the G function of the sources located *INSIDE* the scatterer at the test point;
% % % % % % % %         Mx(i,j) = -1i*omega*(e0*er_out)*scalar_green(in_sources_loc(:,j),test_loc(:,i),n_out,k0,c,OMEGA); % the MINUS sign is due to TE
% % % % % % % %
% % % % % % % %     end
% % % % % % % %     for k = 1:length(out_sources_loc)  % calculate the G function of the sources located *OUTSIDE* the scatterer at the test point, now with with the MINUS sign
% % % % % % % %         Mx(i,j+k) =  1i*omega*(e0*er_in)*scalar_green(out_sources_loc(:,k),test_loc(:,i),n_in,k0,c,OMEGA); % the MINUS sign is due to TE
% % % % % % % %     end
% % % % % % % %     Vex(i) = -H_inc_z(i);
% % % % % % % %
% % % % % % % % end
% % % % % % % %
% % % % % % % %
% % % % % % % % for i=1:length(test_loc) %% calculate nXdeltaH(Z direction)
% % % % % % % %     for j=1:length(in_sources_loc)
% % % % % % % %         [Gx,Gy] = dyiadic_green(in_sources_loc(:,j),test_loc(:,i),n_out,k0);
% % % % % % % %         sol = cross(nhat(i,:),[Gx,Gy,0]);
% % % % % % % %         MHz(i,j) = sol(3);
% % % % % % % %         %         MHz(i,j) = j;
% % % % % % % %     end
% % % % % % % %
% % % % % % % %     for k=1:length(out_sources_loc)
% % % % % % % %         [Gx,Gy] = dyiadic_green(out_sources_loc(:,k),test_loc(:,i),n_in,k0);
% % % % % % % %         sol = cross(nhat(i,:),[Gx,Gy,0]);
% % % % % % % %         MHz(i,j+k) = -1*sol(3);
% % % % % % % %         %         MHz(i,j+k) = -k;
% % % % % % % %     end
% % % % % % % %     sol = cross(nhat(i,:),[E_inc_x(i),E_inc_y(i),0]);
% % % % % % % %     Vhz(i) = -1*sol(3);
% % % % % % % % end


M_total = [MHz;MEt];
% P_inverse_M = ((M_total')*M_total)\(M_total');
% P_inverse_M = (transpose(M_total)*M_total)\(transpose(M_total));

% I_solution = P_inverse_M*(transpose([Vex,Vhz]));
% I_solution = (M_total)\(transpose([Vex,Vhz]));
I_solution = linsolve(M_total,transpose([Vhz,Vet]));
%  I_solution =P_inverse_M*transpose([Vex,Vhz]);



% figure; quiver(test_loc(1,:),test_loc(2,:),direction(1,:),direction(2,:));
% figure; scatter(test_loc(1,:),test_loc(2,:),20,abs(E_inc_z),'filled');
% hold on;
% scatter(in_sources_loc(1,:),in_sources_loc(2,:),20,abs(I_solution(1:length(in_sources_loc))),'filled');
% scatter(out_sources_loc(1,:),out_sources_loc(2,:),20,abs(I_solution((1+length(in_sources_loc)):2*length(in_sources_loc))),'filled');

X = params.X;
Y = params.Y;
loc_line = [X(:),Y(:)]';

H_sol = zeros(1,length(loc_line));


grid_inside_each_sca = ((X).^2 +(Y).^2)<(R*R);
H_sol_final = zeros(length(params.sca_x),length(X(grid_inside_each_sca)));
Ex_sol = H_sol_final;
Ey_sol = H_sol_final;
for tt = 1:length(filaments_total)
    if (filaments_total(tt,4) == 0)
        

        in_loc_line = [X(grid_inside_each_sca)+params.sca_x(filaments_total(tt,3)),Y(grid_inside_each_sca)+params.sca_y(filaments_total(tt,3))]';
        H_sol_final(filaments_total(tt,3),:) = H_sol_final(filaments_total(tt,3),:) + 1i*omega*(e0*er_in)*I_solution(tt)*scalar_green(filaments_total(tt,1:2)',in_loc_line,n_in,k0,c,OMEGA);
        [Gx,Gy] =  dyiadic_green(filaments_total(tt,1:2)',in_loc_line,n_in,k0,c,OMEGA);
        Ex = -Gx*I_solution(tt);
        Ey = -Gy*I_solution(tt);
        Ex_sol(filaments_total(tt,3),:) = Ex_sol(filaments_total(tt,3),:) + Ex;
        Ey_sol(filaments_total(tt,3),:) = Ey_sol(filaments_total(tt,3),:) + Ey;
        


    end
end

tic
for mm = 1 : length(params.sca_x)
    
    in_loc_line = [X(grid_inside_each_sca)+params.sca_x(mm),Y(grid_inside_each_sca)+params.sca_y(mm)]';
    
    for i = 1:length(in_loc_line)
%          E_inc_in_sca(i,:) = double(subs(params.Esym,[xsym,ysym,zsym,k0sym] ,[in_loc_line(1,i),in_loc_line(2,i),0,k0]));
        %             H_inc(i,:) = double(subs(params.Hsym,[xsym,ysym,zsym,k0sym] ,[test_point_total(i,1),test_point_total(i,2),0,k0]));
            H_inc_green_inside(i,:) = scalar_green([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params.n_out,params.k0,params.c,params.OMEGA);

        [E_inc_green_x_temp , E_inc_green_y_temp] = dyiadic_green([params.source_loc_x;params.source_loc_y],[in_loc_line(1,i),in_loc_line(2,i)]',params.n_out,params.k0, c,params.OMEGA);
          E_inc_green_x(i,:) = -1./(1i*params.omega*(params.e0*params.er_out))*E_inc_green_x_temp;
          E_inc_green_y(i,:) = -1./(1i*params.omega*(params.e0*params.er_out))*E_inc_green_y_temp;

    end
%      E_inc_x_inside_sca(mm,:) = E_inc_in_sca(:,1);
%      E_inc_y_inside_sca(mm,:) = E_inc_in_sca(:,2);
    H_inc_z_inside_sca(mm,:) = H_inc_green_inside(:);

    E_inc_x_inside_sca(mm,:) = E_inc_green_x(:);
    E_inc_y_inside_sca(mm,:) = E_inc_green_y(:);
end
toc
% E_inc_x_inside_sca = 0;
% E_inc_y_inside_sca = E_inc_y;


mean_Ex_fil = mean(Ex_sol,2);
mean_Ey_fil = mean(Ey_sol,2);

mean_Ix_fil = mean_Ex_fil*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R);
mean_Iy_fil = mean_Ey_fil*(-1i*params.omega*params.e0*(params.er_in-params.er_out)*pi*R*R);

mean_Ex_inc = mean(E_inc_x_inside_sca,2);
mean_Ey_inc = mean(E_inc_y_inside_sca,2);
mean_Hz_inc = mean(H_inc_z_inside_sca,2);

%           [E_inc_green_x_temp , E_inc_green_y_temp] = dyiadic_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0, c,params.OMEGA);
%           mean_Ex_inc = -1./(1i*params.omega*(params.e0*params.er_out))*E_inc_green_x_temp;
%           mean_Ey_inc = -1./(1i*params.omega*(params.e0*params.er_out))*E_inc_green_y_temp;
%           mean_Hz_inc = scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);






alpha_xx = mean_Ix_fil/mean_Ex_inc;
alpha_yy = mean_Iy_fil/mean_Ey_inc;

alpha_xy = mean_Ix_fil/mean_Ey_inc;
alpha_yx = mean_Iy_fil/mean_Ex_inc;

alpha_mat = [alpha_xx,alpha_xy;alpha_yx,alpha_yy];

output_mat = {mean_Ix_fil,mean_Iy_fil,mean_Ex_inc,mean_Ey_inc,mean_Hz_inc};
% %  alpha_Fil = mean_I_fil/scalar_green([params.source_loc_x;params.source_loc_y],[params.sca_x,params.sca_y]',params.n_out,params.k0,params.c,params.OMEGA);
% alpha_Fil = mean_I_fil/1;


% 
for tt = 1 : length(Hy_sol_final)
   
    final_current_H(tt) = -Hx_sol_final(tt)*in_loc_line(1,tt) -  Hy_sol_final(tt)*in_loc_line(2,tt);
    rotation_current_H(tt) =  -H_inc_x_inside_sca(tt)*in_loc_line(1,tt) -  H_inc_y_inside_sca(tt)*in_loc_line(2,tt);
    Jm_fil_x(tt) = (E_sol_final(tt) - E_inc_z_inside_sca(tt))*in_loc_line(1,tt);
    Jm_fil_y(tt) = (E_sol_final(tt) - E_inc_z_inside_sca(tt))*in_loc_line(2,tt);

end

% figure; quiver(in_loc_line(1,:),in_loc_line(2,:),real(Ex_sol),real(Ey_sol));

% figure; quiver(in_loc_line(1,:),in_loc_line(2,:),imag(Ex_sol-E_inc_x_inside_sca),imag(Ey_sol-E_inc_y_inside_sca));
% figure; quiver(in_loc_line(1,:),in_loc_line(2,:),real(Ex_sol-E_inc_x_inside_sca),real(Ey_sol-E_inc_y_inside_sca));
























%%  syntetic plane wave only
if (false)
    
    X = params.X;
    Y = params.Y;
    loc_line = [X(:),Y(:)]';
    out_loc = (X.*X+Y.*Y)>R*R;
    in_loc = (X.*X+Y.*Y)<=R*R;
    
    out_loc_line = [X(out_loc==1),Y(out_loc==1)]';
    in_loc_line = [X(in_loc==1),Y(in_loc==1)]';
    
    H_sol = zeros(1,length(loc_line));
    Ex_sol = H_sol;
    Ey_sol = H_sol;
    for tt=1:length(in_sources_loc)
        H_sol(out_loc(:)) = H_sol(out_loc(:))+1i*omega*(e0*er_out)*I_solution(tt)*scalar_green(in_sources_loc(:,tt),out_loc_line,n_out,k0,c,OMEGA); % the MINUS sign is due to TE
        [Gx,Gy] =  dyiadic_green(in_sources_loc(:,tt),out_loc_line,n_out,k0,c,OMEGA);
        Ex = -Gx*I_solution(tt);
        Ey = -Gy*I_solution(tt);
        Ex_sol(out_loc(:)) = Ex_sol(out_loc(:)) + Ex;
        Ey_sol(out_loc(:)) = Ey_sol(out_loc(:)) + Ey;
    end
    
    for tt=1:length(out_sources_loc)
        H_sol(in_loc(:)) = H_sol(in_loc(:))+1i*omega*(e0*er_in)*I_solution(tt+length(in_sources_loc))*scalar_green(out_sources_loc(:,tt),in_loc_line,n_in,k0,c,OMEGA);
        [Gx,Gy] =  dyiadic_green(out_sources_loc(:,tt),in_loc_line,n_in,k0,c,OMEGA);
        Ex = -Gx*I_solution(tt+length(in_sources_loc));
        Ey = -Gy*I_solution(tt+length(in_sources_loc));
        Ex_sol(in_loc(:)) = Ex_sol(in_loc(:)) + Ex;
        Ey_sol(in_loc(:)) = Ey_sol(in_loc(:)) + Ey;
        
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
    
    Hsym(3)
    
    
    
    
    
    
    
    
    H_sol_maybe(out_loc(:)) = H_sol(out_loc(:)) +  exp(-1i*k0*out_loc_line(1,:));
    H_sol_maybe(in_loc(:)) = H_sol(in_loc(:));
    
    Ex_sol_maybe(out_loc(:)) = Ex_sol(out_loc(:));
    Ex_sol_maybe(in_loc(:)) = Ex_sol(in_loc(:));
    
    
    
    
    Ey_sol_maybe(out_loc(:)) = Ey_sol(out_loc(:)) -1i/(params.omega*(params.e0*params.er_out))*(-1i*params.k0*exp(-1i*params.k0*out_loc_line(1,:)));
    Ey_sol_maybe(in_loc(:)) = Ey_sol(in_loc(:));
    
    H_sol_sca(out_loc(:)) = H_sol(out_loc(:));
    H_sol_sca(in_loc(:)) = H_sol(in_loc(:))-  exp(-1i*k0*in_loc_line(1,:));
    
    Ex_sol_maybe_sca(out_loc(:)) = Ex_sol(out_loc(:));
    Ex_sol_maybe_sca(in_loc(:)) = Ex_sol(in_loc(:));
    
    Ey_sol_maybe_sca(out_loc(:)) = Ey_sol_maybe(out_loc(:)) -(-1i/(params.omega*(params.e0*params.er_out))*(-1i*params.k0*exp(-1i*params.k0*out_loc_line(1,:))));
    Ey_sol_maybe_sca(in_loc(:)) = Ey_sol_maybe(in_loc(:)) - (-1i/(params.omega*(params.e0*params.er_out))*(-1i*params.k0*exp(-1i*params.k0*in_loc_line(1,:))));
    
    
    %
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol),'filled'); axis equal
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol),'filled');axis equal
    
    
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,real(E_sol_maybe),'filled'); axis equal
    % figure; scatter(loc_line(1,:),loc_line(2,:),20,abs(E_sol_maybe),'filled');axis equal
    
    % figure; imagesc(abs(reshape(H_sol_maybe,length(X),length(X)))); axis equal; title('filaments -- ABS')
    %
    % figure; imagesc(imag(reshape(H_sol_maybe,length(X),length(X)))); axis equal ; title('filaments --IMAG')
    % figure; imagesc(real(reshape(H_sol_maybe,length(X),length(X)))); axis equal ;  title('filaments --REAL')
    H_SOL_st_fil_TE = reshape(H_sol_maybe,length(X),length(X));
    H_sol_sca = reshape(H_sol_sca,length(X),length(X));
    
    Ex_sol_final = reshape(Ex_sol_maybe,length(X),length(X));
    Ey_sol_final = reshape(Ey_sol_maybe,length(X),length(X));
    
    Ex_sol_final_sca = reshape(Ex_sol_maybe_sca,length(X),length(X));
    Ey_sol_final_sca = reshape(Ey_sol_maybe_sca,length(X),length(X));
    
end
end
