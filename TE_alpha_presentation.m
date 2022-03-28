for rho = 1 : size(alpha_te_fil_output,1)
    for OMEGA_r = 1 : size(alpha_te_fil_output,2)
%         temp_alpha = alpha_te_fil_output{rho,OMEGA_r,1};
%         alpha_yy(rho,OMEGA_r) = temp_alpha(2,2);
%         alpha_yx(rho,OMEGA_r) = temp_alpha(1,2);
        field_current{rho,OMEGA_r,1} = field_current_mat_output{rho,OMEGA_r,1};
        
        temp_field_current_1 = cell2mat(field_current{rho,OMEGA_r,1});
        
%         temp_alpha = alpha_te_fil_output{rho,OMEGA_r,2};
%         alpha_xx(rho,OMEGA_r) = temp_alpha(1,1);
%         alpha_xy(rho,OMEGA_r) = temp_alpha(2,1);
        field_current{rho,OMEGA_r,2} = field_current_mat_output{rho,OMEGA_r,2};
        
        temp_field_current_2 = cell2mat(field_current{rho,OMEGA_r,2});
        
       
        field_current{rho,OMEGA_r,3} = field_current_mat_output{rho,OMEGA_r,3};
        
        temp_field_current_3 = cell2mat(field_current{rho,OMEGA_r,3});
        
        Ix_1 = temp_field_current_1(6);
        Iy_1 = temp_field_current_1(7);
        ImZ_1 = temp_field_current_1(1);
        Ex_1 = temp_field_current_1(3);
        Ey_1 = temp_field_current_1(4);
        Hz_1 = temp_field_current_1(5);
        
        Ix_2 = temp_field_current_2(6);
        Iy_2 = temp_field_current_2(7);
        ImZ_2 = temp_field_current_2(1);
        Ex_2 = temp_field_current_2(3);
        Ey_2 = temp_field_current_2(4);
        Hz_2 = temp_field_current_2(5);

        Ix_3 = temp_field_current_3(6);
        Iy_3 = temp_field_current_3(7);
        ImZ_3 = temp_field_current_3(1);
        Ex_3 = temp_field_current_3(3);
        Ey_3 = temp_field_current_3(4);
        Hz_3 = temp_field_current_3(5);
        
        lineq_for_x = [Ex_1 , Ey_1 ,Hz_1; Ex_2 ,Ey_2, Hz_2; Ex_3 ,Ey_3, Hz_3];
        a=linsolve(lineq_for_x,[Ix_1;Ix_2;Ix_3]);
        alpha_exex(rho,OMEGA_r) = a(1);
        alpha_eyex(rho,OMEGA_r) = a(2);
        alpha_hzex(rho,OMEGA_r) = a(3);

        lineq_for_x = [Ex_1 , Ey_1 ,Hz_1; Ex_2 ,Ey_2, Hz_2; Ex_3 ,Ey_3, Hz_3];
        a=linsolve(lineq_for_x,[Iy_1;Iy_2;Iy_3]);
        alpha_exey(rho,OMEGA_r) = a(1);
        alpha_eyey(rho,OMEGA_r) = a(2);
        alpha_hzey(rho,OMEGA_r) = a(3);
        
        
        lineq_for_x = [Ex_1 , Ey_1 ,Hz_1; Ex_2 ,Ey_2, Hz_2; Ex_3 ,Ey_3, Hz_3];
        a=linsolve(lineq_for_x,[ImZ_1;ImZ_2;ImZ_3]);
        alpha_exhz(rho,OMEGA_r) = a(1);
        alpha_eyhz(rho,OMEGA_r) = a(2);
        alpha_hzhz(rho,OMEGA_r) = a(3);
        
        
        
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



% figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_yy),'x-','LineWidth',2); title ('\alpha_{\theta\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_yx),'x-','LineWidth',2); title ('\alpha_{\rho\theta}$$','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1);plot(params.OMEGA_vec/params.omega,abs(alpha_xx),'x-','LineWidth',2); title ('\alpha_{\rho\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs(alpha_xy),'x-','LineWidth',2); title ('\alpha_{\theta\rho}$$','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% legend (num2str(params.shift_vec/params.lambda'),'LineWidth',2);
% 
% 
% figure;
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyy),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxy),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1); plot(params.OMEGA_vec/params.omega,abs(alpha_lineqxx),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyx),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% 
% figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyy),'*-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,real(alpha_lineqxy),'*-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1);plot(params.OMEGA_vec/params.omega,real(alpha_lineqxx),'*-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,real(alpha_lineqyx),'*-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% 
% figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyy),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxy),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex'); xlabel ('\Omega/\omega'); grid on; grid minor;
% subplot(2,2,1);plot(params.OMEGA_vec/params.omega,imag(alpha_lineqxx),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag(alpha_lineqyx),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex');xlabel ('\Omega/\omega'); grid on; grid minor;
% 
% 
% 
% 
% figure; subplot(2,2,4); (surf(abs(alpha_lineqyy))); title ('E at $$\hat{\theta} , \alpha_{\theta\theta}$$','Interpreter','Latex'); grid on; grid minor;
% subplot(2,2,2); (surf(abs(alpha_lineqxy))); title ('E at $$\hat{\theta} , \alpha_{\rho\theta}$$','Interpreter','Latex'); grid on; grid minor;
% subplot(2,2,1);(surf(abs(alpha_lineqxx))); title ('E at $$\hat{\rho} , \alpha_{\rho\rho}$$','Interpreter','Latex');grid on; grid minor;
% subplot(2,2,3); (surf(abs(alpha_lineqyx))); title ('E at $$\hat{\rho} , \alpha_{\theta\rho}$$','Interpreter','Latex');grid on; grid minor;
% 

% % % figure;
% % % subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_eyey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_eyex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_exex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % %  lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% % %  title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,3,4); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_exey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hzex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hzey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % figure;
% % % subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_eyey),'x-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_eyex),'x-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_exex),'x-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % %  lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% % %  title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,3,4); p4 =plot(params.OMEGA_vec/params.omega,real(alpha_exey),'x-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hzex),'x-','LineWidth',2); title ('real($$\alpha_{\rho H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hzey),'x-','LineWidth',2); title ('real($$\alpha_{\theta H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % 
% % % figure;
% % % subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_eyey),'x-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_eyex),'x-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_exex),'x-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % %  lgd1 = legend (num2str(params.shift_vec'/params.lambda));
% % %  title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,3,4); p4 =plot(params.OMEGA_vec/params.omega,imag(alpha_exey),'x-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_hzex),'x-','LineWidth',2); title ('imag($$\alpha_{\rho H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,3,6); plot(params.OMEGA_vec/params.omega,imag(alpha_hzey),'x-','LineWidth',2); title ('imag($$\alpha_{\theta H}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % 
% % % 
% % % figure;
% % % subplot(2,2,4); plot(params.OMEGA_vec/params.omega,abs(alpha_eyey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,abs(alpha_eyex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,1); plot(params.OMEGA_vec/params.omega,abs(alpha_exex),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % % % title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_exey),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,real(alpha_eyey),'*-','LineWidth',2); title ('real($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,real(alpha_eyex),'*-','LineWidth',2); title ('real($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega'); grid on; grid minor;
% % % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,real(alpha_exex),'*-','LineWidth',2); title ('real($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % % % title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,real(alpha_exey),'*-','LineWidth',2); title ('real($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % 
% % % figure; subplot(2,2,4); plot(params.OMEGA_vec/params.omega,imag(alpha_eyey),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,2); plot(params.OMEGA_vec/params.omega,imag(alpha_eyex),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % subplot(2,2,1);plot(params.OMEGA_vec/params.omega,imag(alpha_exex),'s-','LineWidth',2); title ('imag($$\alpha_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% % % % lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% % % % title(lgd1,'\rho_c [\lambda]')
% % % subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag(alpha_exey),'s-','LineWidth',2); title ('imag($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

% 
% figure;
% subplot(2,2,1); plot(params.OMEGA_vec/params.omega,real((alpha_lineqxx + alpha_lineqyx)./alpha_lineqxx(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,2); plot(params.OMEGA_vec/params.omega,atan2(real(alpha_lineqxx),imag(alpha_lineqyx)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,imag((alpha_lineqxx + alpha_lineqyx)./alpha_lineqxx(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,atan2(imag(alpha_lineqxx),imag(alpha_lineqyx)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% 
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs((alpha_lineqxx + alpha_lineqyx)./alpha_lineqxx(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,angle(alpha_lineqxx + alpha_lineqyx)-angle(alpha_lineqxx(3,6)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% 
% subplot(2,2,3); plot(params.OMEGA_vec/params.omega,abs((alpha_lineqyy + alpha_lineqxy)./alpha_lineqyy(6)-1),'x-','LineWidth',2); title ('abs($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% lgd1 = legend (num2str(params.shift_vec'/params.lambda),'LineWidth',2);
% subplot(2,2,4); plot(params.OMEGA_vec/params.omega,angle(alpha_lineqyy + alpha_lineqxy)-angle(alpha_lineqyy(3,6)),'x-','LineWidth',2); title ('\angle($$\alpha_{\rho\rho} + \alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% title(lgd1,'\rho_c [\lambda]')
% subplot(2,2,3); p4 =plot(params.OMEGA_vec/params.omega,abs(alpha_lineqyx),'x-','LineWidth',2); title ('abs($$\alpha_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% 
% 
% 

% 
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_hzhz),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_eyhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_exhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_hzex),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_hzey),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%  
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  
%  title(lgd1,'\rho_c [\lambda]')
%  
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_hzhz),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_eyhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_exhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_hzex),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_hzey),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
%  
% 
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%  
%  figure;
% subplot(2,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_hzhz),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_eyhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_exhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_hzex),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% subplot(2,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_hzey),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
% 
% 
% lgd1 = legend (num2str(params.shift_vec'/params.lambda));
%  title(lgd1,'\rho_c [\lambda]')
%  
%  
 
 
 
 
 figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,abs(alpha_hzhz),'x-','LineWidth',2); title ('abs($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,abs(alpha_eyhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,abs(alpha_exhz),'x-','LineWidth',2); title ('abs($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,abs(alpha_hzex),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,abs(alpha_hzey),'x-','LineWidth',2); title ('abs($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,abs(alpha_exex),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,abs(alpha_exey),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,abs(alpha_eyex),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,abs(alpha_eyey),'x-','LineWidth',2); title ('abs($$\alpha^{ee}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));
 
 title(lgd1,'\rho_c [\lambda]')
  
 figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,real(alpha_hzhz),'d-','LineWidth',2); title ('real($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,real(alpha_eyhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,real(alpha_exhz),'d-','LineWidth',2); title ('real($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,real(alpha_hzex),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,real(alpha_hzey),'d-','LineWidth',2); title ('real($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,real(alpha_exex),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,real(alpha_exey),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,real(alpha_eyex),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,real(alpha_eyey),'d-','LineWidth',2); title ('real($$\alpha^{ee}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));
 
 title(lgd1,'\rho_c [\lambda]')
 
 figure;
subplot(3,3,9); plot(params.OMEGA_vec/params.omega,imag(alpha_hzhz),'+-','LineWidth',2); title ('imag($$\alpha^{mm}_{zz}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,8); plot(params.OMEGA_vec/params.omega,imag(alpha_eyhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,7); plot(params.OMEGA_vec/params.omega,imag(alpha_exhz),'+-','LineWidth',2); title ('imag($$\alpha^{me}_{z\rho}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,3); plot(params.OMEGA_vec/params.omega,imag(alpha_hzex),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\rho z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,6); plot(params.OMEGA_vec/params.omega,imag(alpha_hzey),'+-','LineWidth',2); title ('imag($$\alpha^{em}_{\theta z}$$)','Interpreter','Latex', 'FontSize', 16);xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,1); plot(params.OMEGA_vec/params.omega,imag(alpha_exex),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\rho\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,2); plot(params.OMEGA_vec/params.omega,imag(alpha_exey),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\rho\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,4); plot(params.OMEGA_vec/params.omega,imag(alpha_eyex),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\theta\rho}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;
subplot(3,3,5); plot(params.OMEGA_vec/params.omega,imag(alpha_eyey),'+-','LineWidth',2); title ('imag($$\alpha^{ee}_{\theta\theta}$$)','Interpreter','Latex', 'FontSize', 16); xlabel ('\Omega/\omega', 'FontSize', 16); grid on; grid minor;

lgd1 = legend (num2str(params.shift_vec'/params.lambda));
 
 title(lgd1,'\rho_c [\lambda]')
