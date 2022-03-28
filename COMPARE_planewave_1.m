%% Compare All Methods -
% 3 scatterer
% radius [1/100, 1/50] wavelength
% located {[0,0], [0.2,0] ,[-0.1,0.1]}
% plane wave
% mu = 1
% epsilon = 4 , 11
% TM Only
% No rotation




eps_vec = [4,7,10];
radius_vec = [1/100, 1/50];



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

% params.radius = 0.01*1e-6; % [meters]

params.source_loc_x = lambda*[-1];
params.source_loc_y = lambda*[0.5];

params.sca_x = lambda*[ 0 ];
params.sca_y = lambda*[ 0 ];

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

    %    Parameters of Filamets;
    params.R_out=1.2; %params.radius*param.R_out
    params.R_in=0.8;
    
    params.N_filaments=100; %on each side
    params.N_testpoints=200;%must be >=param.N_filaments

    params.max_m = 50;
%% START SIMULATION
    % Hitting Field
    E_inc_z = exp(-1i*k0*params.X);   
    E_inc_x = 0*E_inc_z;
    E_inc_y = 0*E_inc_z;
    
    H_inc_y = 1i/(omega*(mu0*mr_out))*(-1i*k0*exp(-1i*k0*params.X));
    H_inc_x = 0*(H_inc_y);
    H_inc_z = 0*(H_inc_y);
    
    %symbolic
    syms xsym ysym zsym k0sym
    Esym = [0, 0,exp(-1i*k0sym*xsym)]; % plane wave
    Xsym = [xsym ysym zsym];
    Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);
    % end

for i = 1:length(eps_vec)
    params.er_in = eps_vec(i);
    params.n_in = sqrt(params.mr_in*params.er_in);
    params.n_out = sqrt(params.mr_out*params.er_out);
    
    for j = 1:length(radius_vec)
        params.radius = radius_vec(j)*lambda;
        
        %MIE SERIES
        tic
        EVAL_ALL_RESULTS;
              
    end
end
now_str = datestr(now,'mmmm_dd_yyyy_HH_MM_SS');
save (['Planewave_1Sca_results_',now_str]);

PRINT_RESULTS;

        

