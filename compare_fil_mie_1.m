%% COMPARE TM & TE with filaments and Mie Series
close all;

declare_params;

te = 1;
tm = 1;
toplot = 1;

disp(params);

%% if TE, calculate the Incident Hz, Et fields
if (te)
    % set the Ez_incident and the Ht_incident
    H_inc_z = exp(-1i*k0*params.X);
    H_inc_x = 0*H_inc_z;
    H_inc_y = 0*H_inc_z;
    
    E_inc_y = -1i/(omega*(e0*er_out))*(-1i*k0*exp(-1i*k0*params.X));
    E_inc_x = 0*(E_inc_y);
    E_inc_z = 0*(E_inc_y);
    
    [H_SOL_st_mie,H_sca_Mie] = Mie_Series_TE(params,H_inc_z);
    
    syms xsym ysym zsym k0sym
    %     Esym = [0, 0,exp(-1i*k0sym*xsym)];
    %     Xsym = [xsym ysym zsym];
    %     Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);
    
    Hsym = [0, 0,exp(-1i*k0sym*xsym)];
    Xsym = [xsym ysym zsym];
    Esym = 1i/(omega*(e0*er_out))*curl(Hsym,Xsym);

    
    [H_SOL_st_filaments,H_sca_fil] = filaments_TE(params,H_inc_z,E_inc_x,E_inc_y,Hsym);
    if (toplot)
        figure;
        subplot(1,3,1); imagesc(params.x,params.y,abs(H_SOL_st_filaments)); title('filaments TE - ABS H'); axis equal; colorbar;
        subplot(1,3,2); imagesc(params.x,params.y,abs(H_SOL_st_mie)); title('Mie TE - ABS H'); axis equal; colorbar;
        subplot(1,3,3); imagesc(params.x,params.y,abs(H_SOL_st_mie)-abs(H_SOL_st_filaments)); title('Diff TE H'); axis equal; colorbar;
        
        figure;
        subplot(1,3,1); imagesc(params.x,params.y,real(H_SOL_st_filaments)); title('filaments TE - real H'); axis equal; colorbar;
        subplot(1,3,2); imagesc(params.x,params.y,real(H_SOL_st_mie)); title('Mie TE - real H'); axis equal; colorbar;
        subplot(1,3,3); imagesc(params.x,params.y,real(H_SOL_st_mie)-real(H_SOL_st_filaments)); title('Diff TE real H'); axis equal; colorbar;
        
        
        figure; subplot(1,3,1); imagesc(params.x,params.y,abs(H_sca_fil));
        title('filaments TE, Only Scattered field ABS'); axis equal; colorbar; grid minor;
        subplot(1,3,2); imagesc(params.x,params.y,real(H_sca_fil));
        title('filaments TE, Only Scattered field REAL'); axis equal; colorbar; grid minor;
        subplot(1,3,3); imagesc(params.x,params.y,imag(H_sca_fil));
        title('filaments TE, Only Scattered field IMAG'); axis equal; colorbar; grid minor;
        
        
        figure; subplot(1,3,1); imagesc(params.x,params.y,abs(H_sca_Mie));
        title('Mie TE, Only Scattered field ABS'); axis equal; colorbar; grid minor;
        subplot(1,3,2); imagesc(params.x,params.y,real(H_sca_Mie));
        title('Mie TE, Only Scattered field REAL'); axis equal; colorbar; grid minor;
        subplot(1,3,3); imagesc(params.x,params.y,imag(H_sca_Mie));
        title('Mie TE, Only Scattered field IMAG'); axis equal; colorbar; grid minor;
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
    
    [E_SOL_st_mie,E_sca_Mie] = Mie_Series_TM(params,E_inc_z);
    
    syms xsym ysym zsym k0sym
    Esym = [0, 0,exp(-1i*k0sym*xsym)];
    Xsym = [xsym ysym zsym];
    Hsym = -1i/(omega*(mu0*mr_out))*curl(Esym,Xsym);
    
    
    [E_SOL_st_filaments, E_sca_fil] = filaments_TM(params,E_inc_z,H_inc_x,H_inc_y,Esym);
    
    if (toplot)
        
        figure;
        subplot(1,3,1); imagesc(params.x,params.y,abs(E_SOL_st_filaments)); title('filaments TM - ABS E'); axis equal; colorbar;
        subplot(1,3,2); imagesc(params.x,params.y,abs(E_SOL_st_mie)); title('Mie TM - ABS E'); axis equal; colorbar;
        subplot(1,3,3); imagesc(params.x,params.y,abs(E_SOL_st_mie)-abs(E_SOL_st_filaments)); title('Diff TM E'); axis equal; colorbar;
        
        figure;
        subplot(1,3,1); imagesc(params.x,params.y,real(E_SOL_st_filaments)); title('filaments TM - real E'); axis equal; colorbar;
        subplot(1,3,2); imagesc(params.x,params.y,real(E_SOL_st_mie)); title('Mie TM - real E'); axis equal; colorbar;
        subplot(1,3,3); imagesc(params.x,params.y,real(E_SOL_st_mie)-real(E_SOL_st_filaments)); title('Diff TM real E'); axis equal; colorbar;
        
        
        figure; subplot(1,3,1); imagesc(params.x,params.y,abs(E_sca_fil));
        title('filaments TM, Only Scattered field ABS'); axis equal; colorbar; grid minor;
        subplot(1,3,2); imagesc(params.x,params.y,real(E_sca_fil));
        title('filaments TM, Only Scattered field REAL'); axis equal; colorbar; grid minor;
        subplot(1,3,3); imagesc(params.x,params.y,imag(E_sca_fil));
        title('filaments TM, Only Scattered field IMAG'); axis equal; colorbar; grid minor;
        
        
        figure; subplot(1,3,1); imagesc(params.x,params.y,abs(E_sca_Mie));
        title('Mie TM, Only Scattered field ABS'); axis equal; colorbar; grid minor;
        subplot(1,3,2); imagesc(params.x,params.y,real(E_sca_Mie));
        title('Mie TM, Only Scattered field REAL'); axis equal; colorbar; grid minor;
        subplot(1,3,3); imagesc(params.x,params.y,imag(E_sca_Mie));
        title('Mie TM, Only Scattered field IMAG'); axis equal; colorbar; grid minor;
        
    end
    
end

disp(params)

figure;  imagesc(params.x,params.y,abs(E_SOL_st_filaments-E_SOL_st_mie)./(abs(E_SOL_st_filaments)+abs(E_SOL_st_mie)));
title('filaments TM, ABS E_{fil} - ABS E_{mie} / Sum(ABS)'); axis equal; colorbar; grid minor;

figure;  imagesc(params.x,params.y,abs(H_SOL_st_filaments-H_SOL_st_mie)./(abs(H_SOL_st_filaments)+abs(H_SOL_st_mie)));
title('filaments TE, ABS H_{fil} - ABS H_{mie} / Sum(ABS)'); axis equal; colorbar; grid minor;

figure; subplot(1,3,1); imagesc(params.x,params.y,abs(E_sca_fil));
title('filaments TM, Only Scattered field ABS'); axis equal; colorbar; grid minor;
subplot(1,3,2); imagesc(params.x,params.y,real(E_sca_fil));
title('filaments TM, Only Scattered field REAL'); axis equal; colorbar; grid minor;
subplot(1,3,3); imagesc(params.x,params.y,imag(E_sca_fil));
title('filaments TM, Only Scattered field IMAG'); axis equal; colorbar; grid minor;