for rho = 1 : size(E_SOL_st_alpha_MoM,1)
    for OMEGA_r = 1 : size(E_SOL_st_alpha_MoM,2)

    alpha_1(rho,OMEGA_r) =  E_SOL_st_alpha_MoM{rho,OMEGA_r,1};
    alpha_2(rho,OMEGA_r) =  E_SOL_st_alpha_MoM{rho,OMEGA_r,2};
    alpha_3(rho,OMEGA_r) =  E_SOL_st_alpha_MoM{rho,OMEGA_r,3};
    
    E_mom_1 (rho,OMEGA_r) = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,1};
    E_mom_2 (rho,OMEGA_r) = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,2};
    E_mom_3 (rho,OMEGA_r) = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,3};

   OMEGA(rho,OMEGA_r) = params.OMEGA_vec(OMEGA_r)/params.omega;
   rho_c(rho,OMEGA_r) = params.shift_vec(rho);
   
        field_current{rho,OMEGA_r,1} = field_current_mat_output{rho,OMEGA_r,1};
        temp_field_current_1 = cell2mat(field_current{rho,OMEGA_r,1});
        
         Iz_1 = temp_field_current_1(1);
         Iez_1 = temp_field_current_1(2);
        Hx_1 = temp_field_current_1(3);
        Hy_1 = temp_field_current_1(4);
        Ez_1 = temp_field_current_1(5);
        Im_x_1 = temp_field_current_1(6);
        Im_y_1 = temp_field_current_1(7);
         Ez_inc_1 = temp_field_current_1(8);

        
        field_current{rho,OMEGA_r,2} = field_current_mat_output{rho,OMEGA_r,2};
        temp_field_current_2 = cell2mat(field_current{rho,OMEGA_r,2});
        
        Iz_2 = temp_field_current_2(1);
        Iez_2 = temp_field_current_2(2);
        Hx_2 = temp_field_current_2(3);
        Hy_2 = temp_field_current_2(4);
        Ez_2 = temp_field_current_2(5);
        Im_x_2 = temp_field_current_2(6);
        Im_y_2 = temp_field_current_2(7);
         Ez_inc_2 = temp_field_current_2(8);

       
        field_current{rho,OMEGA_r,3} = field_current_mat_output{rho,OMEGA_r,3};
        temp_field_current_3 = cell2mat(field_current{rho,OMEGA_r,3});
        
        Iz_3 = temp_field_current_3(1);
        Iez_3 = temp_field_current_3(2);
        Hx_3 = temp_field_current_3(3);
        Hy_3 = temp_field_current_3(4);
        Ez_3 = temp_field_current_3(5);
        Im_x_3 = temp_field_current_3(6);
        Im_y_3 = temp_field_current_3(7);
            Ez_inc_3 = temp_field_current_3(8);

        
    alpha_1_fil(rho,OMEGA_r) =  Iez_1/Ez_1;
    alpha_2_fil(rho,OMEGA_r) =  Iez_2/Ez_2;
    alpha_3_fil(rho,OMEGA_r) =  Iez_3/Ez_3;
   
    E_1_fil(rho,OMEGA_r) =  Ez_inc_1;
    E_2_fil(rho,OMEGA_r) =  Ez_inc_2;
    E_3_fil(rho,OMEGA_r) =  Ez_inc_3;
      
end

end

figure;
subplot(1,3,1);
surf(OMEGA,rho_c,abs(alpha_1)); 
title ('|\alpha_{TM} 1| of a 1/100 \lambda radius cylinder, \epsilon_r = 11.4 - MoM');
xlabel ('\Omega/\omega');
ylabel ('\rho_c [\lambda]'); 

subplot(1,3,2);
surf(OMEGA,rho_c,abs(alpha_2)); 
title ('|\alpha_{TM} 2 | of a 1/100 \lambda radius cylinder, \epsilon_r = 11.4 - MoM');
xlabel ('\Omega/\omega');
ylabel ('\rho_c [\lambda]'); 


subplot(1,3,3);
surf(OMEGA,rho_c,abs(alpha_3)); 
title ('|\alpha_{TM} 3 | of a 1/100 \lambda radius cylinder, \epsilon_r = 11.4 - MoM');
xlabel ('\Omega/\omega');
ylabel ('\rho_c [\lambda]'); 


figure; plot(abs(alpha_1),'.-');
hold on;
plot(abs(alpha_1_fil),'*-');
figure; plot(abs(alpha_2),'*-'); 
hold on;
plot(abs(alpha_2_fil),'*-');
figure; plot(abs(alpha_3),'d-');
hold on;
plot(abs(alpha_3_fil),'*-');



figure; plot(abs(E_mom_1),'.-');
hold on;
plot(abs(E_1_fil),'*-');
figure; plot(abs(E_mom_2),'*-'); 
hold on;
plot(abs(E_2_fil),'*-');
figure; plot(abs(E_mom_3),'d-');
hold on;
plot(abs(E_3_fil),'*-');



figure; subplot(2,1,1);plot(abs(alpha_1)./abs(alpha_1(1,:)),'.-');
hold on;plot(abs(alpha_1_fil)./abs(alpha_1_fil(1,:)),'*--');

figure; subplot(2,1,1);plot(abs(alpha_2)./abs(alpha_2(1,:)),'.-');
hold on;plot(abs(alpha_2_fil)./abs(alpha_2_fil(1,:)),'*--');

figure; subplot(2,1,1);plot(abs(alpha_3)./abs(alpha_3(1,:)),'.-');
subplot(2,1,2); plot(abs(alpha_3_fil)./abs(alpha_3_fil(1,:)),'*--');


figure; subplot(2,1,1);plot(abs(E_mom_1)./abs(E_mom_1(1,1)),'.-');
hold on;plot(abs(E_1_fil)./abs(E_1_fil(1,1)),'s--');

figure; subplot(2,1,1);plot(abs(E_mom_2)./abs(E_mom_2(1,1)),'.-');
hold on;plot(abs(E_2_fil)./abs(E_2_fil(1,1)),'d--');

figure; subplot(2,1,1);plot(abs(E_mom_3)./abs(E_mom_3(1,1)),'.-');
subplot(2,1,2); plot(abs(E_3_fil)./abs(E_3_fil(1,1)),'*--');


