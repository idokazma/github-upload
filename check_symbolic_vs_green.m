% params.H_inc_x = reshape(params.H_inc_x,length(params.X),[]);
% params.H_inc_y = reshape(params.H_inc_y,length(params.X),[]);

H_inc_x_SYM = reshape(H_inc(:,1),length(params.X),[]);
H_inc_y_SYM = reshape(H_inc(:,2),length(params.X),[]);

E_inc_z_SYM = reshape(E_inc(:,3),length(params.X),[]);
figure; imagesc(abs(params.E_inc_z-E_inc_z_SYM));colorbar;


figure;
subplot(1,2,1); imagesc(abs(H_inc_x_SYM)); colorbar;set(gca,'YDir','normal')
subplot(1,2,2); imagesc(abs(params.H_inc_x));colorbar;set(gca,'YDir','normal')
figure;imagesc(abs(params.H_inc_x-H_inc_x_SYM));colorbar;set(gca,'YDir','normal')
figure;imagesc(abs(params.H_inc_y-H_inc_y_SYM));colorbar;set(gca,'YDir','normal')

figure;imagesc(imag(params.H_inc_y));colorbar;set(gca,'YDir','normal')
figure;imagesc(imag(H_inc_y_SYM));colorbar;set(gca,'YDir','normal')

figure;imagesc(imag(params.H_inc_y-H_inc_y_SYM));colorbar;set(gca,'YDir','normal')
figure;imagesc(imag(params.H_inc_x-H_inc_x_SYM));colorbar;set(gca,'YDir','normal')

figure;imagesc(real(params.H_inc_y-H_inc_y_SYM));colorbar;set(gca,'YDir','normal')
figure;imagesc(real(params.H_inc_x-H_inc_x_SYM));colorbar;set(gca,'YDir','normal')


figure;imagesc(abs(params.H_inc_y-H_inc_y_SYM));colorbar;set(gca,'YDir','normal')



figure;
subplot(1,2,1); imagesc(abs(H_inc_y_SYM));colorbar;
subplot(1,2,2); imagesc(abs(params.H_inc_y));colorbar;

figure;
subplot(1,2,1); imagesc(real(H_inc_x_SYM));colorbar;
subplot(1,2,2); imagesc(real(params.H_inc_x));colorbar;

figure;
subplot(1,2,1); imagesc(imag(H_inc_x_SYM));colorbar;
subplot(1,2,2); imagesc(imag(params.H_inc_x));colorbar;

figure;
subplot(1,2,1); imagesc(real(H_inc_y_SYM));colorbar;
subplot(1,2,2); imagesc(real(params.H_inc_y));colorbar;