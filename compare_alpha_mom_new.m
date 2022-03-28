rotations = (-1:0.2:1)*1e-3;
% rotations = (0:5:50)*1e-7;

Mom = [];
Pol = []
for i=1:length(rotations)

      Mom(:,i) = MoM_tri_1_6(-rotations(i),VogelArrayXY)
tic
     Pol(:,i) = RotatingArray_2D_TM(VogelArrayXY',rotations(i));
toc
rotations(i)
close all;
end

figure; plot(rotations,abs(Pol)','-s');
hold on; plot(rotations,abs(Mom)','-o');
grid on; grid minor
title ('Polarization and MoM Polarization Currents for 5 scatteres');
xlabel ('\Omega/\omega');
ylabel ('|Polarization Current|');

figure; plot(rotations, abs(Mom)'./abs(Pol)')

figure; plot(rotations, abs(Mom)'-abs(Pol)')

figure; plot(abs(Pol(1,:))' -abs(Mom(1,:))' ,'-.');


figure; plot(rotations,abs(Pol)');
for ii = 1:length(rotations)
    resu(ii,:) = sum((abs(Pol)./(abs(Pol(:,101))) - abs(Pol(:,ii))./(abs(Pol(:,101)))).^2,1);
    resu_2(ii,:) = sum((abs(Pol)' - abs(Pol(:,ii))').^2,2);

    resu_coh(ii,:) = sum(abs(((Pol)' - (Pol(:,ii))')).^2,2);

end
figure; imagesc(rotations,rotations,(resu)); colormap('jet');
figure; imagesc(rotations,rotations,(resu_2)); colormap('jet');

figure; imagesc(rotations,rotations,db(resu_coh)); colormap('jet');


figure; imagesc(rotations,1:2000,abs(Pol)./(abs(Pol(:,76)))); colormap('jet');
figure; imagesc(rotations,1:2000,abs(Pol)./(abs(Pol(:,76)))); colormap('jet');

Ratios = abs(Pol)./(abs(Pol(:,26)));
ratios_diff = diff(Ratios');
figure; imagesc(rotations(1:end-1),1:2000,ratios_diff'); colormap('jet');

hold on; plot(rotations,abs(Mom)');
grid on; grid minor
title ('Polarization and MoM Polarization Currents for 5 scatteres');
xlabel ('\Omega/\omega');
ylabel ('|Polarization Current|');