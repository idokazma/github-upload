
for i = 1:size(PARAMETERS_PLANE,1)
    for j = 1:size(PARAMETERS_PLANE,2)

        params = PARAMETERS_PLANE{i,j};
        E_SOL_st_MoM = E_SOL_st_MoM_PLANE{i,j};
        E_SOL_st_filaments_mul = E_SOL_st_filaments_mul_PLANE{i,j};
        
        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,abs(E_SOL_st_filaments_mul-params.E_inc_z)); title(['FIL Sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,abs(E_SOL_st_MoM-params.E_inc_z)); title(['MoM sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,abs((E_SOL_st_MoM)-(E_SOL_st_filaments_mul))); title(['DIFF, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;

        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,real(E_SOL_st_filaments_mul-params.E_inc_z)); title(['FIL Sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,real(E_SOL_st_MoM-params.E_inc_z)); title(['MoM sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,real((E_SOL_st_MoM)-(E_SOL_st_filaments_mul))); title(['DIFF, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        
        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,imag(E_SOL_st_filaments_mul-params.E_inc_z)); title(['FIL Sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,imag(E_SOL_st_MoM-params.E_inc_z)); title(['MoM sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,imag((E_SOL_st_MoM)-(E_SOL_st_filaments_mul))); title(['DIFF, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;

        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,real(E_SOL_st_filaments_mul-params.E_inc_z)); title(['FIL Sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,real(E_SOL_st_MoM-params.E_inc_z)); title(['MoM sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,real((E_SOL_st_MoM)-(E_SOL_st_filaments_mul))); title(['DIFF, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        
        
        
        figure;
        subplot(1,3,1); imagesc(params.x, params.y ,imag(E_SOL_st_filaments_mul)); title(['FIL Sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,2); imagesc(params.x, params.y ,imag(E_SOL_st_MoM)); title(['MoM sca, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;
        subplot(1,3,3); imagesc(params.x, params.y ,imag((E_SOL_st_MoM)-(E_SOL_st_filaments_mul))); title(['DIFF, r=' num2str(params.radius/params.lambda),'\lambda epsR=',num2str(params.er_in)]); axis equal; colorbar;draw_sca_line;

        
    end
end