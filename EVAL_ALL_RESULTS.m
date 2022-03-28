%% EVAL_ALL_RESULTS

if (tm)

[E_SOL_st_MoM, mean_E_MoM,mean_I_MoM,effective_radius] = MoM_tri_1_6(params);
% mean_E_MoM
% mean_I_MoM
        fprintf ('done MoM\n')
close all;
params.radius = effective_radius;

[Pvec, E_SOL_st_POL] = RotatingArray_2D_TM(params, params.E_inc_z, params.Esym);
% E_SOL_st_POL
% Pvec
   fprintf ('done POL\n')
end


if (params.is_plane_wave == 1 && length(params.sca_x)== 1 && (params.sca_x)== 0 && (params.sca_y)== 0)
[E_SOL_st_mie, mean_E_mie] = Mie_Series_TM(params,params.E_inc_z);
mean_E_mie = mean(mean_E_mie)
        fprintf('done MIE\n')
end
[E_SOL_st_filaments_mul, mean_E_fil] = filaments_TM_multiple(params,params.Esym);
mean_E_fil
         fprintf ('done FIL\n')


     
        
        
        
  PARAMETERS_PLANE{i,j} = params;
try
    E_SOL_st_mie_PLANE{i,j} = E_SOL_st_mie;
    E_SOL_st_mie_MEAN_PLANE{i,j} = mean_E_mie;
catch
    E_SOL_st_mie_PLANE{i,j} = -1; 
end

E_SOL_st_filaments_mul_PLANE{i,j} = E_SOL_st_filaments_mul;
E_SOL_st_filaments_mul_I_PLANE{i,j} = mean_E_fil;

E_SOL_st_MoM_PLANE{i,j} = E_SOL_st_MoM;
E_SOL_st_MoM_MEAN_PLANE{i,j} = mean_E_MoM;
E_SOL_st_MoM_I_MEAN_PLANE{i,j} = mean_I_MoM;

E_SOL_st_POL_PLANE{i,j} = E_SOL_st_POL;
E_SOL_st_POL_I_PLANE{i,j} = Pvec;
toc

fprintf ('~~~~ DONE eps %d, radius %d ~~~~ \n', i/length(eps_vec)*100, j/length(radius_vec)*100)