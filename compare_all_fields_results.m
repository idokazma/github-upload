for d = 1 : (size(E_SOL_st_POL_PLANE,1))

for mm = 1 : (size(E_SOL_st_POL_PLANE,2))
    
%    Fil_E(mm,:) = cell2mat(E_SOL_st_filaments_mul_E_PLANE(d,mm));
%    MoM_E(mm,:) = cell2mat(E_SOL_st_MoM_MEAN_PLANE(d,mm)).';
%    Pol_E(mm,:) = cell2mat(E_SOL_st_POL_PLANE(d,mm));
%    
%    Fil_I(mm,:) = cell2mat(E_SOL_st_filaments_mul_I_PLANE(d,mm));
%    MoM_I(mm,:) = cell2mat(E_SOL_st_MoM_I_MEAN_PLANE(d,mm)).';
%    Pol_I(mm,:) = cell2mat(E_SOL_st_POL_I_PLANE(d,mm));

   Fil_alpha(mm,d) = cell2mat(E_SOL_st_filaments_mul_alpha(d,mm,1));
   MoM_alpha(mm,d) = cell2mat(E_SOL_st_alpha_MoM(d,mm,1)).';
   Pol_alpha(mm,d) = cell2mat(E_SOL_st_alpha_POL(d,mm,1));
   OMEGA(mm,d) = params.OMEGA_vec(mm)/params.omega;
   rho_c(mm,d) = params.shift_vec(d);
   
      Fil_alpha_plane(mm,d) = cell2mat(E_SOL_st_filaments_mul_alpha(d,mm,2));
   MoM_alpha_plane(mm,d) = cell2mat(E_SOL_st_alpha_MoM(d,mm,2)).';
   Pol_alpha_plane(mm,d) = cell2mat(E_SOL_st_alpha_POL(d,mm,2));
   OMEGA(mm,d) = params.OMEGA_vec(mm)/params.omega;
   rho_c(mm,d) = params.shift_vec(d);
   
end

end
figure; surf(OMEGA,rho_c,abs(Fil_alpha)); 
title ('|\alpha_{TM}| of a 1/100 \lambda radius cylinder, \epsilon_r = 11.4 - FIL');
xlabel ('\Omega/\omega');
ylabel ('\rho_c [\lambda]');

figure; surf(OMEGA,rho_c,abs(MoM_alpha)); 
title ('|\alpha_{TM}| of a 1/100 \lambda radius cylinder, \epsilon_r = 11.4 - MoM');
xlabel ('\Omega/\omega');
ylabel ('\rho_c [\lambda]');

% figure;   plot(abs(Fil_E)./abs(MoM_E),'.-b'); hold on; plot(abs(MoM_E)./abs(MoM_E),'.-r'); plot(abs(Pol_E)./abs(MoM_E),'.-g');
% figure;   plot(abs(Fil_E),'.-b'); hold on; plot(abs(MoM_E),'.-r'); plot(abs(Pol_E),'.-g');
% figure;   plot(abs(Fil_E),'.-b'); hold on; plot(abs(Pol_E),'.-g');
% 
% figure;   plot(abs(Fil_E./(Fil_E(1))-1),'.-b'); hold on;
% plot(abs(MoM_E./(MoM_E(1))-1),'.-r'); 
% plot(abs(Pol_E./(Pol_E(1))-1),'.-g'); 
% 
% 
% figure;   plot(abs(Fil_E./(Fil_E(1))),'.-b'); hold on;
% plot(abs(MoM_E./(MoM_E(1))),'.-r'); 
% plot(abs(Pol_E./(Pol_E(1))),'.-g');
% 

figure;   plot(params.OMEGA_vec/params.omega,abs(Fil_alpha),'.-'); hold on;
   plot(params.OMEGA_vec/params.omega,abs(Fil_alpha_plane),'+-'); hold on;

figure;   plot(params.OMEGA_vec/params.omega,abs(MoM_alpha),'.-'); hold on;
   plot(params.OMEGA_vec/params.omega,abs(MoM_alpha_plane),'+-'); hold on;
%    
% figure;plot(params.OMEGA_vec/params.omega,abs(MoM_alpha),'.-r'); 
% plot(params.OMEGA_vec/params.omega,abs(Pol_alpha),'.-g');
% title([num2str(params.N_filaments),' filaments on each side']); 
% legend('Fil','MoM','Mie Static');
% 
% figure;  plot(params.OMEGA_vec/params.omega,abs(Fil_alpha./(Fil_alpha(1))),'.-b'); hold on;
% plot(params.OMEGA_vec/params.omega,abs(MoM_alpha./(MoM_alpha(1))),'.-r'); 
% plot(params.OMEGA_vec/params.omega,abs(Pol_alpha./(Pol_alpha(1))),'.-g'); 
% title([num2str(params.N_filaments),' filaments on each side']); 
% legend('Fil','MoM','Mie Static');
% 
% figure;   plot(params.OMEGA_vec/params.omega,abs(Fil_alpha./(Fil_alpha(1))-1),'.-b'); hold on;
% plot(params.OMEGA_vec/params.omega,abs(MoM_alpha./(MoM_alpha(1))-1),'.-r'); 
% plot(params.OMEGA_vec/params.omega,abs(Pol_alpha./(Pol_alpha(1))-1),'.-g');
% title([num2str(params.N_filaments),' filaments on each side']); 
% legend('Fil','MoM','Mie Static');
% 
% 
% figure; plot(rho_c,abs(Fil_alpha),'.-'); hold on; 
% plot(rho_c,abs(MoM_alpha),'.-');
% plot(rho_c,abs(Pol_alpha),'.-');
% 
% figure; plot(abs(Fil_alpha),'.-'); hold on; 
% plot(abs(MoM_alpha),'.-');
% plot(abs(Pol_alpha),'.-');
% 
% ratio_fil = abs(Fil_alpha./(Fil_alpha(1))-1);
% ratio_mom = abs(MoM_alpha./(MoM_alpha(1))-1);
% OMEGA_vec = params.OMEGA_vec;
% rho_c = params.lambda*((1:20) + 0);
% 
% figure; plot(rho_c,ratio_fil,'.-')
% hold on; plot(rho_c,ratio_mom,'.-r')
% 
% 
% figure; plot(ratio_fil,'.-')
% hold on; plot(ratio_mom,'.-r')
% 
% 
% ratio_fil_pol = abs(Fil_alpha./(Pol_alpha(1))-1);
% ratio_mom_pol = abs(MoM_alpha./(Pol_alpha(1))-1);
% OMEGA_vec = params.OMEGA_vec;
% rho_c = params.lambda*((1:20) );
% 
% figure; plot(rho_c,ratio_fil_pol,'.-')
% hold on; plot(rho_c,ratio_mom_pol,'.-r')
% % hold on; plot(abs(MoM_E)/abs(MoM_E(1)),'.-g');
% % plot(abs(Pol_E)/abs(Pol_E(1)),'.-r');
% % 
% % figure;   plot(abs(Fil_I),'.-b'); hold on; plot(abs(MoM_I),'.-r'); plot(abs(Pol_I),'.-g');
% % 
% % figure; plot(abs(MoM_I)/abs(MoM_I(1)),'.-g');
% 
% % transpose(cell2mat(E_SOL_st_filaments_mul_E_PLANE))
% % transpose(cell2mat(E_SOL_st_MoM_MEAN_PLANE))
% % transpose(cell2mat(E_SOL_st_POL_PLANE))
% % 
% % 
% % 
% % 
% % 
% % transpose(cell2mat(E_SOL_st_POL_I_PLANE))
% % transpose(cell2mat(E_SOL_st_MoM_I_MEAN_PLANE))
