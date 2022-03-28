%% COMPARE TM & TE with filaments and Mie Series
te = 0;
tm = 1;
toplot = 1;

start_time = datestr(now);
disp (start_time);
%declare GRID to calculate the fields in;

params.len_n = 601;
params.wid_n = 601;
params.len = 0.51e-6;
params.wid = 0.51e-6;

params.y = linspace(0,params.len,params.len_n)- params.len/2;
params.x = linspace(0,params.wid,params.wid_n)- params.wid/2;

[params.X,params.Y] = meshgrid(params.x,params.y);


% params.len_n_mom = 201;
% params.wid_n_mom = 201;
% params.len_mom = 2.1*1e-6;
% params.wid_mom = 2.1*1e-6;

params.y_mom = linspace(0,params.len,params.len_n*1)- params.len/2;
params.x_mom = linspace(0,params.wid,params.wid_n*1)- params.wid/2;




%% State the physical parameters fo the problem (epsilon(r), mu(r), radius...)
params.lambda = 1e-6; %[meters]
lambda = params.lambda;

params.radius = 0.01*1e-6; % [meters]

params.sca_x = lambda*[ 0];
params.sca_y = lambda*[ 0 ];

params.source_loc_x = lambda*[-1];
params.source_loc_y = lambda*[0.5];

params.is_plane_wave = 1;

params.e0= 8.85418781762039*1e-12;
e0 = params.e0;

params.mu0 = 4*pi*1e-7;
mu0 = params.mu0;

params.mr_in = 1;
mr_in = params.mr_in;
params.mr_out = 1;
mr_out = params.mr_out;
params.er_in = 11;
er_in = params.er_in;
params.er_out = 1;
er_out = params.er_out;

params.c=1/sqrt(params.e0*params.mu0);
c =params.c;

params.f = params.c/params.lambda;
f = params.f;

params.k0 = 2*pi/params.lambda;
k0 = params.k0;

params.omega = params.f*2*pi;
omega = params.omega;

params.n_in = sqrt(params.mr_in*params.er_in);
params.n_out = sqrt(params.mr_out*params.er_out);

params.OMEGA = 0*1e-3*omega;

%% Choose if TE or TM

%% if TE, calculate the Incident Hz, Et fields
if (te)
    % set the Ez_incident and the Ht_incident
    H_inc_z = exp(-1i*k0*params.X);
    H_inc_x = 0*H_inc_z;
    H_inc_y = 0*H_inc_z;
    
    E_inc_y = -1i/(omega*(e0*er_out))*(-1i*k0*exp(-1i*k0*params.X));
    E_inc_x = 0*(E_inc_y);
    E_inc_z = 0*(E_inc_y);
    
    params.max_m = 50;
    [H_SOL_st_mie] = Mie_Series_TE(params,H_inc_z);
    
    syms xsym ysym zsym k0sym
    %     Esym = [0, 0,exp(-1i*k0sym*xsym)];
    %     Xsym = [xsym ysym zsym];
    %     Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);
    
    Hsym = [0, 0,exp(-1i*k0sym*xsym)];
    Xsym = [xsym ysym zsym];
    Esym = 1i/(omega*(e0*er_out))*curl(Hsym,Xsym);

    %    Parameters of Filamets;
    params.R_out=1.2; %params.radius*param.R_out
    params.R_in=0.8;
    
    params.N_filaments=100; %on each side
    params.N_testpoints=300;%must be >=param.N_filaments
    
    [H_SOL_st_filaments] = filaments_TE(params,H_inc_z,E_inc_x,E_inc_y,Hsym);
    if (toplot)
        figure;
        subplot(1,3,1); imagesc(params.x,params.y,abs(H_SOL_st_filaments)); title('filaments TE - ABS H'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x,params.y,abs(H_SOL_st_mie)); title('Mie TE - ABS H'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x,params.y,abs(H_SOL_st_mie)-abs(H_SOL_st_filaments)); title('Diff TE H'); axis equal; colorbar;draw_sca_line;
        
        figure;
        subplot(1,3,1); imagesc(params.x,params.y,real(H_SOL_st_filaments)); title('filaments TE - real H'); axis equal; colorbar;
        subplot(1,3,2); imagesc(params.x,params.y,real(H_SOL_st_mie)); title('Mie TE - real H'); axis equal; colorbar;
        subplot(1,3,3); imagesc(params.x,params.y,real(H_SOL_st_mie)-real(H_SOL_st_filaments)); title('Diff TE real H'); axis equal; colorbar;
    end
end

%% if TM, calculate the Incident Ez, Ht fields
if (tm)
    % set the Ez_incident and the Ht_incident
    E_inc_z = exp(-1i*k0*params.X);   
    E_inc_x = 0*E_inc_z;
    E_inc_y = 0*E_inc_z;
    
    H_inc_y = 1i/(omega*(mu0*mr_out))*(-1i*k0*exp(-1i*k0*params.X));
    H_inc_x = 0*(H_inc_y);
    H_inc_z = 0*(H_inc_y);
    %%%
    params.max_m = 50;
      [E_SOL_st_mie] = Mie_Series_TM(params,E_inc_z);
    fprintf ('done MIE')
    syms xsym ysym zsym k0sym
    Esym = [0, 0,exp(-1i*k0sym*ysym)]; % plane wave
%     ez = 1i/4 * besselh (0 , 1 , k0sym *sqrt(xsym.^2+ysym.^2));
%     Esym = [0, 0, ez];

    Xsym = [xsym ysym zsym];
    Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);

    %    Parameters of Filamets;
    params.R_out=1.2; %params.radius*param.R_out
    params.R_in=0.8;
    
    params.N_filaments=200; %on each side
    params.N_testpoints=300;%must be >=param.N_filaments
    
%     [E_SOL_st_filaments] = filaments_TM(params,E_inc_z,H_inc_x,H_inc_y,Esym);

     [E_SOL_st_filaments_mul] = filaments_TM_multiple(params,E_inc_z,H_inc_x,H_inc_y,Esym);
         fprintf ('done FIL')

    [E_SOL_st_MoM] = MoM_tri_1_6(params);
        fprintf ('done MoM')
     [E_SOL_st_POL] = 0;    
        

    
    if (toplot)
%       figure;
%       subplot(1,3,1); imagesc(abs(E_SOL_st_filaments)); title('filaments TM - ABS E'); axis equal; colorbar;
%       subplot(1,3,2); imagesc(abs(E_SOL_st_mie)); title('Mie TM - ABS E'); axis equal; colorbar;
%       subplot(1,3,3); imagesc(abs(E_SOL_st_mie)-abs(E_SOL_st_filaments)); title('Diff TM E'); axis equal; colorbar;
        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,real(E_SOL_st_filaments_mul-E_inc_z)); title('filaments Sca'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,real(E_SOL_st_MoM-E_inc_z)); title('MoM sca'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,real((E_SOL_st_MoM-E_inc_z)-(E_SOL_st_filaments_mul-E_inc_z))); title('Diff'); axis equal; colorbar;draw_sca_line;

        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,abs(E_SOL_st_filaments_mul-E_inc_z)); title('filaments Sca'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,abs(E_SOL_st_MoM-E_inc_z)); title('MoM sca'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,abs((E_SOL_st_MoM-E_inc_z)./(E_SOL_st_filaments_mul-E_inc_z))); title('Ratio'); axis equal; colorbar;draw_sca_line;
        
        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,abs(E_SOL_st_filaments_mul)); title('filaments TM ABS E'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,abs(E_SOL_st_mie)); title('Mie ABS E'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,abs(E_SOL_st_MoM)); title('MOM ABS E'); axis equal; colorbar;draw_sca_line;
        
        figure;
        subplot(1,2,1); imagesc(params.x, params.y ,abs(E_SOL_st_filaments_mul-E_SOL_st_mie)); title('filaments TM - Mie TM ABS E'); axis equal; colorbar;draw_sca_line;
        subplot(1,2,2); imagesc(params.x, params.y ,abs(E_SOL_st_MoM-E_SOL_st_mie)); title('MOM - Mie ABS E'); axis equal; colorbar;draw_sca_line;

        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,real(E_SOL_st_filaments_mul)); title('filaments TM ABS E'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,real(E_SOL_st_mie)); title('Mie ABS E'); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,real(E_SOL_st_MoM)); title('MOM ABS E'); axis equal; colorbar;draw_sca_line;
        
        figure;
        subplot(1,2,1); imagesc(params.x, params.y ,real(E_SOL_st_filaments_mul-E_SOL_st_mie)); title('filaments TM - Mie TM real E'); axis equal; colorbar;draw_sca_line;
        subplot(1,2,2); imagesc(params.x, params.y ,real(E_SOL_st_MoM-E_SOL_st_mie)); title('MOM - Mie real E'); axis equal; colorbar;draw_sca_line;
        
        
        figure;
        imagesc(params.x, params.y ,abs(E_SOL_st_filaments_mul-E_SOL_st_MoM)); title('filaments TM - Mie TM ABS E'); axis equal; colorbar;draw_sca_line;
        
%         figure;
%         subplot(1,3,1); imagesc(real(E_SOL_st_filaments)); title('filaments TM - real E'); axis equal; colorbar;
%         subplot(1,3,2); imagesc(real(E_SOL_st_mie)); title('Mie TM - real E'); axis equal; colorbar;
%         subplot(1,3,3); imagesc(real(E_SOL_st_mie)-real(E_SOL_st_filaments)); title('Diff TM real E'); axis equal; colorbar;
    end
    
end

disp(params)

fprintf ('~~~Start time was: %s \n',start_time);
fprintf ('~~~~~~Finished at: %s\n', datestr(now));