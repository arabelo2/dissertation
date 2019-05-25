loadingDataOfRectangularArrayPiston
cd E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\
processa
figure(50)
[rowsOfMaxes colsOfMaxes] = find(Ppp == max(max(Ppp)));
hold on; plot(lambda*(1:length(Ppp(rowsOfMaxes, :))), abs(Ppp(rowsOfMaxes, :))./max(max(Ppp)))
%%%% hold on; plot(lambda*(1:length(Ppp(round(length(x(1, :))/ 2), :))), Ppp(round(length(x(1, :))/ 2), :)./max(max(Ppp)), 'b')
picoapicointerpolado = resample(picoapico', 8, 1)';
hold on; plot(0.5*(1:length(picoapico(:, 60))), picoapico(:, 60)/max(picoapico(:, 60)))
maxPosition = round(find(picoapicointerpolado == max(max(picoapicointerpolado)))/size(picoapicointerpolado, 1))+1;
hold on; plot(0.5*(1:length(picoapicointerpolado(:, maxPosition))), picoapicointerpolado(:, maxPosition)/max(max(picoapicointerpolado)))
title ('Pressão acústica ao longo do eixo | fc = 5MHz e fs = 100 MHz | F = 40 mm | DADOS0H80.DAT')
grid on
grid minor
legend('Simulated pulse', 'Experimental pulse')
ylim([0,1])
legend('Simulated pulse DADOS0H80', 'Experimental pulse')
loadingDataOfRectangularArrayPiston
plot(lambda*(1:length(Ppp(round(length(x(1, :))/ 2), :))), Ppp(round(length(x(1, :))/ 2), :)./max(max(Ppp)))
legend('Simulated pulse DADOS0H80', 'Experimental pulse', 'Simulated pulse DADOS0H0' )
title ('Pressão acústica ao longo do eixo | fc = 5MHz e fs = 100 MHz | F = 40 mm')
xlabel('z (mm)')
close all
