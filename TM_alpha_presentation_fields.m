%fields only

for rho = 1 : size(field_current_mat_output,1)
    for OMEGA_r = 1 : size(field_current_mat_output,2)
        
        field_current{rho,OMEGA_r,1} = field_current_mat_output{rho,OMEGA_r,1};
        temp_field_current_1 = cell2mat(field_current{rho,OMEGA_r,1});
        
        field_current{rho,OMEGA_r,2} = field_current_mat_output{rho,OMEGA_r,2};
        temp_field_current_2 = cell2mat(field_current{rho,OMEGA_r,2});
        
        field_current{rho,OMEGA_r,3} = field_current_mat_output{rho,OMEGA_r,3};
        temp_field_current_3 = cell2mat(field_current{rho,OMEGA_r,3});
        
        
        Ez_sol_MoM_1 = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,1};
        Ez_sol_MoM_2 = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,2};
        Ez_sol_MoM_3 = E_SOL_st_MoM_MEAN_PLANE{rho,OMEGA_r,3};
        
        Ez_inc_MoM_1 = E_SOL_st_MoM_MEAN_PLANE_E0{rho,OMEGA_r,1};
        Ez_inc_MoM_2 = E_SOL_st_MoM_MEAN_PLANE_E0{rho,OMEGA_r,2};
        Ez_inc_MoM_3 = E_SOL_st_MoM_MEAN_PLANE_E0{rho,OMEGA_r,3};
        
        
        Iz_1 = temp_field_current_1(1);
        Iez_1 = temp_field_current_1(2);
        Hx_1 = temp_field_current_1(3);
        Hy_1 = temp_field_current_1(4);
        Ez_1 = temp_field_current_1(5);
        Im_x_1 = temp_field_current_1(6);
        Im_y_1 = temp_field_current_1(7);
        Ez_sol_1 = temp_field_current_1(8);
        Hx_sol_1 = temp_field_current_1(9);
        Hy_sol_1 = temp_field_current_1(10);

        
        Iz_2 = temp_field_current_2(1);
        Iez_2 = temp_field_current_2(2);
        Hx_2 = temp_field_current_2(3);
        Hy_2 = temp_field_current_2(4);
        Ez_2 = temp_field_current_2(5);
        Im_x_2 = temp_field_current_2(6);
        Im_y_2 = temp_field_current_2(7);
        Ez_sol_2 = temp_field_current_2(8);
        Hx_sol_2 = temp_field_current_2(9);
        Hy_sol_2 = temp_field_current_2(10);
        
        Iz_3 = temp_field_current_3(1);
        Iez_3 = temp_field_current_3(2);
        Hx_3 = temp_field_current_3(3);
        Hy_3 = temp_field_current_3(4);
        Ez_3 = temp_field_current_3(5);
        Im_x_3 = temp_field_current_3(6);
        Im_y_3 = temp_field_current_3(7);
        Ez_sol_3 = temp_field_current_3(8);
        Hx_sol_3 = temp_field_current_3(9);
        Hy_sol_3 = temp_field_current_3(10);
        
        log_ez(rho,OMEGA_r,1) = Ez_1;
        log_ez(rho,OMEGA_r,2) = Ez_2;
        log_ez(rho,OMEGA_r,3) = Ez_3;
        
        log_alpha(rho,OMEGA_r,1) = Iez_1/Ez_1;
        log_alpha(rho,OMEGA_r,2) = Iez_2/Ez_2;
        log_alpha(rho,OMEGA_r,3) = Iez_3/Ez_3;

        log_I(rho,OMEGA_r,1) = Iez_1;
        log_I(rho,OMEGA_r,2) = Im_x_1;
        log_I(rho,OMEGA_r,3) = Im_y_1;
        
        lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
        a=linsolve((lineq_for_x),[Ez_sol_1;Ez_sol_2;Ez_sol_3]);
        alpha_ezhx(rho,OMEGA_r) = a(1);
        alpha_ezhy(rho,OMEGA_r) = a(2);
        alpha_ezez(rho,OMEGA_r) = a(3);
        
%         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=inv(lineq_for_x)*([Iz_1;Iz_2;Iz_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
        
%         lineq_for_x = [0 , 0 ,Ez_1; 0 ,0, Ez_2; 0 ,0, Ez_3];
%         a=linsolve(lineq_for_x,[Iz_1;Iz_2;Iz_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
     
%         
%         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=linsolve(lineq_for_x,[Iez_1;Iez_2;Iez_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
%         
%             lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=linsolve(lineq_for_x,[Iz_1-Iez_1;Iz_2-Iez_2;Iz_3-Iez_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
        %
        %         lineq_for_x = [Ez_1; Ez_2;Ez_3];
        a=linsolve(lineq_for_x,[Hx_sol_1;Hx_sol_2;Hx_sol_3]);
        alpha_hxhx(rho,OMEGA_r) = a(1);
        alpha_hxhy(rho,OMEGA_r) = a(2);
        alpha_hxez(rho,OMEGA_r) = a(3);
        
        %         lineq_for_x = [Ez_1; Ez_2;Ez_3];
        a=linsolve(lineq_for_x,[Hy_sol_1;Hy_sol_2;Hy_sol_3]);
        alpha_hyhx(rho,OMEGA_r) = a(1);
        alpha_hyhy(rho,OMEGA_r) = a(2);
        alpha_hyez(rho,OMEGA_r) = a(3);
        
%         %% fields alpha
%         a=linsolve(lineq_for_x,[Hx_sol_1;Hx_sol_2;Hx_sol_3]);
%         alpha_hxhx(rho,OMEGA_r) = a(1);
%         alpha_hxhy(rho,OMEGA_r) = a(2);
%         alpha_hxez(rho,OMEGA_r) = a(3);
%         
%         %         lineq_for_x = [Ez_1; Ez_2;Ez_3];
%         a=linsolve(lineq_for_x,[Hy_sol_1;Hy_sol_2;Hy_sol_3]);
%         alpha_hyhx(rho,OMEGA_r) = a(1);
%         alpha_hyhy(rho,OMEGA_r) = a(2);
%         alpha_hyez(rho,OMEGA_r) = a(3);
%         
%         
%           lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
%         a=linsolve(lineq_for_x,[Ez_sol_1;Ez_sol_2;Ez_sol_3]);
%         alpha_ezhx(rho,OMEGA_r) = a(1);
%         alpha_ezhy(rho,OMEGA_r) = a(2);
%         alpha_ezez(rho,OMEGA_r) = a(3);
%         
        
        
        
        
        
        
        
        
        
        
        
        lineq_for_x = [Ez_1;  Ez_2;  Ez_3];
        alpha_TM(rho,OMEGA_r)=linsolve(lineq_for_x,[Ez_sol_1;Ez_sol_2;Ez_sol_3]);
        
        %                 lineq_for_x = [Ez_1;  Ez_2;  Ez_3];
        %         alpha_TM(rho,OMEGA_r)=linsolve(lineq_for_x,[Iz_1;Iz_2;Iz_3]);
        %
        lineq_for_x = [Ez_inc_MoM_1;  Ez_inc_MoM_2;  Ez_inc_MoM_3];
        alpha_TM_MoM(rho,OMEGA_r)=linsolve(lineq_for_x,[Ez_sol_MoM_1;Ez_sol_MoM_2;Ez_sol_MoM_3]);
        
        
        %         alpha_TM_MoM(rho,OMEGA_r) = E_SOL_st_alpha_MoM{rho,OMEGA_r,1};

        %         alpha_hxez(rho,OMEGA_r) = a(1);
        %         alpha_hyez(rho,OMEGA_r) = a(2);
        %         alpha_ezez(rho,OMEGA_r) = a(3);
        
        
        %          b(rho,OMEGA_r) = mean([Iez_1 / Ez_1,Iez_2 / Ez_2,Iez_3 / Ez_3]);
        %          lineq_for_x = [Hx_1 , Hy_1 ; Hx_2 ,Hy_2; Hx_3 ,Hy_3];
        %         a=linsolve(lineq_for_x,[Iz_1-Iez_1;Iz_2-Iez_2;Iz_3-Iez_3]);
        %         alpha_lineqhx(rho,OMEGA_r) = a(1);
        %         alpha_lineqhy(rho,OMEGA_r) = a(2);
        %         alpha_lineqez(rho,OMEGA_r) = b(rho,OMEGA_r);
        
        %         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
        %         a=inv(lineq_for_x)*[Iz_1;Iz_2;Iz_3];
        %         alpha_lineqhx(rho,OMEGA_r) = a(1);
        %         alpha_lineqhy(rho,OMEGA_r) = a(2);
        %         alpha_lineqez(rho,OMEGA_r) = a(3);
        %
        
        %         lineq_for_x = [Hx_1 , Hy_1 ,Ez_1; Hx_2 ,Hy_2, Ez_2; Hx_3 ,Hy_3, Ez_3];
        %         a=linsolve(lineq_for_x,[Iy_1;Iy_2;Iy_3]);
        %         alpha_lineqyx(rho,OMEGA_r) = a(1);
        %         alpha_lineqyy(rho,OMEGA_r) = a(2);
        %         alpha_lineqyh(rho,OMEGA_r) = a(3);
        
        %         lineq_for_x = [Ex_1 , Ey_1; Ex_2 ,Ey_2; Ex_3 ,Ey_3];
        %         a=linsolve(lineq_for_x,[Iy_1;Iy_2;Iy_3]);
        %         alpha_lineqyx(rho,OMEGA_r) = a(1);
        %         alpha_lineqyy(rho,OMEGA_r) = a(2);
        %         alpha_lineqyh(rho,OMEGA_r) = 0;
        %
        %         lineq_for_x = [Ex_1 , Ey_1; Ex_2 ,Ey_2; Ex_3 ,Ey_3];
        %         a=linsolve(lineq_for_x,[Ix_1;Ix_2;Ix_3]);
        %         alpha_lineqxx(rho,OMEGA_r) = a(1);
        %         alpha_lineqxy(rho,OMEGA_r) = a(2);
        %         alpha_lineqxh(rho,OMEGA_r) = 0;
        
        %         lineq_for_y = [Ex_1 , Ey_1; Ex_2 ,Ey_2];
        %         a=linsolve(lineq_for_y,[Iy_1;Iy_2]);
        %         alpha_lineqyx(rho,OMEGA_r) = a(1);
        %         alpha_lineqyy(rho,OMEGA_r) = a(2);
        
        
    end
end


%
% % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_yy),'x-','LineWidth',2); title ('\alpha_{\theta\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_yx),'x-','LineWidth',2); title ('\alpha_{\rho\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,abs(alpha_xx),'x-','LineWidth',2); title ('\alpha_{\rho\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs(alpha_xy),'x-','LineWidth',2); title ('\alpha_{\theta\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % legend (num2str(params.shift_vec/params.lambda'),'LineWidth',2);
% %
% %
% % figure;
% % subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyy),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxy),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxx),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % title(lgd1,'\rho_c [\lambda]')
% % subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyx),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% %
% % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyy),'*-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqxy),'*-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,real(alpha_lineqxx),'*-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % title(lgd1,'\rho_c [\lambda]')
% % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyx),'*-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% %
% % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyy),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxy),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxx),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % title(lgd1,'\rho_c [\lambda]')
% % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyx),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% %
% %
% %
% %
% % figure; subplot(2,2,4); (surf(abs(alpha_lineqyy))); title ('E at $$\hat{\theta} , \alpha_{\theta\theta}$$','Interpreter','Latex'); grid on; grid minor;
% % subplot(2,2,2); (surf(abs(alpha_lineqxy))); title ('E at $$\hat{\theta} , \alpha_{\rho\theta}$$','Interpreter','Latex'); grid on; grid minor;
% % subplot(2,2,1);(surf(abs(alpha_lineqxx))); title ('E at $$\hat{\rho} , \alpha_{\rho\rho}$$','Interpreter','Latex');grid on; grid minor;
% % subplot(2,2,3); (surf(abs(alpha_lineqyx))); title ('E at $$\hat{\rho} , \alpha_{\theta\rho}$$','Interpreter','Latex');grid on; grid minor;
% %
%
% figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez),'x-','LineWidth',2); title ('abs($$\alpha_{Ez}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhy),'x-','LineWidth',2); title ('abs($$\alpha_{H\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhx),'x-','LineWidth',2); title ('abs($$\alpha_{H\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Ix),'x-','LineWidth',2); title ('abs($$\alpha_{mEz\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Iy),'x-','LineWidth',2); title ('abs($$\alpha_{mEz\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez),'x-','LineWidth',2); title ('real($$\alpha_{Ez}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhy),'x-','LineWidth',2); title ('real($$\alpha_{H\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhx),'x-','LineWidth',2); title ('real($$\alpha_{H\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Ix),'x-','LineWidth',2); title ('real($$\alpha_{mEz\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Iy),'x-','LineWidth',2); title ('real($$\alpha_{mEz\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez),'x-','LineWidth',2); title ('imag($$\alpha_{Ez}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhy),'x-','LineWidth',2); title ('imag($$\alpha_{H\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhx),'x-','LineWidth',2); title ('imag($$\alpha_{H\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Ix),'x-','LineWidth',2); title ('imag($$\alpha_{mEz\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Iy),'x-','LineWidth',2); title ('imag($$\alpha_{mEz\theta}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%
%
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhy),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqhx),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Ix),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqez_Iy),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhy),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqhx),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Ix),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_lineqez_Iy),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhy),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqhx),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Ix),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqez_Iy),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%
%  figure;
% subplot(2,3,1); plot((alpha_lineqez),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,(alpha_lineqhy),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,(alpha_lineqhx),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,(alpha_lineqez_Ix),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,(alpha_lineqez_Iy),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%
%
%  figure;
%  plot(alpha_lineqez,'o-'); grid on; grid minor;
%   figure;
%  plot(alpha_lineqhx,'.-'); grid on; grid minor;
%
%
%  figure;
% subplot(1,3,1); plot(params.OMEGA_vec/params.omega,abs(b),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%

figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_ezez),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhy),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_ezhx),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hxez),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hyez),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhx),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_hxhy),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhx),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hyhy),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,real(alpha_ezez),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,real(alpha_ezhy),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,real(alpha_ezhx),'d-','LineWidth',2); title ('real($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hxez),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hyez),'d-','LineWidth',2); title ('real($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_hxhx),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_hxhy),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_hyhx),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_hyhy),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')

figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,imag(alpha_ezez),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,imag(alpha_ezhy),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,imag(alpha_ezhx),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_hxez),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,imag(alpha_hyez),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_hxhx),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_hxhy),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_hyhx),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_hyhy),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));

title(lgd1,'\rho_c [\lambda]')


figure;
plot(params.OMEGA_vec/params.omega, abs(alpha_TM),'x-'); hold on;
figure;
plot(params.OMEGA_vec/params.omega, abs(alpha_TM_MoM),'d-');
grid on;

