% cd C:\Users\rabeloal\Documents\MATLAB\TA
% loadingDataOfRectangularArrayPiston
% Ppp = load('fs40.mat'); Ppp = Ppp.Ppp;
[rowsOfMaxes colsOfMaxes] = find(Ppp == max(max(Ppp)));
xMin = 2;
xMax = 68;
yMin = -25;
yMax = 25;

% cd E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\
cd C:\Users\rabeloal\Documents\MATLAB\Experience\var0
processa
picoapicointerpolado = resample(picoapico', 10, 1)';
x_maxPosition = round(find(picoapicointerpolado == max(max(picoapicointerpolado)))/size(picoapicointerpolado, 1))+1;
x_axisCentralPosition = ceil(size(picoapicointerpolado, 2)/2);

figure(50)
% Experimental signal - Pure
hold on; plot(linspace(xMin, xMax, size(picoapico, 1)), picoapico(:, 84)/max(max(picoapico)), 'ko')

% Experimental signal - Interpolated
%hold on; plot(linspace(xMin, xMax, size(picoapicointerpolado, 1)), picoapicointerpolado(:,x_axisCentralPosition)/max(max(picoapicointerpolado)))
hold on; plot(linspace(xMin, xMax, size(picoapicointerpolado, 1)), picoapicointerpolado(:,x_maxPosition)/max(max(picoapicointerpolado)), 'b.')
% Simulated signal
hold on; plot(z(length(find(z < xMin*1e-3))+1:end)*1e3, abs(Ppp(rowsOfMaxes, length(find(z < xMin*1e-3))+1:end))./max(max(Ppp)), 'r')
grid on
grid minor
ylim([0,1])
xlim([0,105])
legend('Sinal experimental', 'Sinal experimental interpolado na direção x', 'Sinal teórico')
title (['Pressão acústica ao longo do eixo para F = 40 mm com deflexão {\theta_0 = }' num2str(THETA), '^{\circ}'])
xlabel('z (mm)', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel(' Amplitude Normalizado', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
set(gca,'FontSize',25);

%legend('Experimental pulse', 'Simulated pulse' )

figure(60)
[peak, y_maxPosition ]= max(picoapicointerpolado(:, x_maxPosition));
hold on; plot(linspace(yMin, yMax, size(picoapico, 2)), picoapico(y_maxPosition, :)/max(max(picoapico)), 'ko')
hold on; plot(linspace(yMin, yMax, size(picoapicointerpolado, 2)), picoapicointerpolado(y_maxPosition,:)/max(max(picoapicointerpolado)), 'b.')
hold on; plot(x*1e3, abs(Ppp(:, length(find(z < (F*cos(deg2rad(THETA))))) + 1)./max(max(Ppp))), 'r')
%hold on; plot(x*1e3, abs(Ppp(:, colsOfMaxes))./max(max(Ppp)), 'ro')
grid on
grid minor
ylim([0,1])
xlim([-31,31])
%title ('Perfil da pressão acústica na posição z = 80 mm')
title (['Perfil da pressão acústica para F = 40 mm com deflexão {\theta_0 = }' num2str(THETA), '^{\circ}'])
legend('Sinal experimental', 'Sinal experimental interpolado na direção x', 'Sinal teórico')
xlabel('x (mm)', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('Amplitude Normalizada', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
set(gca,'FontSize',25);

figure(61)
[peak, y_maxPosition ]= max(picoapicointerpolado(:, x_maxPosition));

xe_angle = rad2deg(asin(1e-3*linspace(yMin, yMax, size(picoapico, 2))/F));
xei_angle = rad2deg(asin(1e-3*linspace(yMin, yMax, size(picoapicointerpolado, 2))/F));
xt_angle = rad2deg(asin(x/F));

ye_db = 20*log10(picoapico(y_maxPosition, :)/max(max(picoapico)));
yei_db = 20*log10(picoapicointerpolado(y_maxPosition,:)/max(max(picoapicointerpolado)));
yt_db = 20*log10(abs(Ppp(:, length(find(z < (F*cos(deg2rad(THETA))))) + 1)./max(max(Ppp))));

hold on; plot(xe_angle, ye_db, 'ko')
hold on; plot(xei_angle, yei_db, 'b.')
hold on; plot(xt_angle, yt_db, 'r')
xlim([-60,60])
grid on
grid minor
%title ('Perfil da pressão acústica na posição z = 80 mm')
title (['Perfil da pressão acústica para F = 40 mm com deflexão {\theta_0 = }' num2str(THETA), '^{\circ}'])
legend('Sinal experimental', 'Sinal experimental interpolado na direção x', 'Sinal teórico')
xlabel('$${\theta} ({^\circ})$$', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('dB = 20log10(Ppp)', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
set(gca,'FontSize',25);

figure(65)
Ht = Ppp(:, length(find(z < (F*cos(deg2rad(THETA))))) + 1)./max(max(Ppp))
hold on ; plot(10*log10(Ht))
He = picoapico(y_maxPosition,:)/max(max(picoapico))
hold on; plot(linspace(1, length(Ht), length(He)), 10*log10(He))
He = picoapicointerpolado(y_maxPosition,:)/max(max(picoapicointerpolado))
hold on; plot(linspace(1, length(Ht), length(He)), 10*log10(He))
grid on
grid minor

figure(50)
hold on; plot(linspace(xMin, xMax, size(picoapicointerpolado, 1)), picoapico(:, 60)/max(max(picoapico)), 'b')
hold on; plot(z(length(find(z < xMin*1e-3))+1:end)*1e3, abs(Ppp(rowsOfMaxes, length(find(z < xMin*1e-3))+1:end))./max(max(Ppp)), 'r')
hold on; plot(linspace(xMin, xMax, size(picoapicointerpolado, 1)), picoapico(:, 60)/max(max(picoapico)), 'k.')
grid on
grid minor
ylim([0,1])
title ('Foco = 40 mm e x = 0 mm')
legend('Sinal experimental (puro)', 'Sinal teórico')
xlabel('z (mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel(' Amplitude Normalizado', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
set(gca,'FontSize',20);


figure(70)
hold on; plot(0.5*(1:length(picoapico(:, 60))), picoapico(:, 60)/max(picoapico(:, 60)))
hold on; plot(0.5*(1:length(picoapicointerpolado(:, x_maxPosition))), picoapicointerpolado(:, x_maxPosition)/max(max(picoapicointerpolado)))
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])

figure(80)
hold on; plot(0.5*(1:length(picoapico(:, 50))), picoapico(:, 50)/max(max(picoapico)))
hold on; plot(0.5*(1:length(picoapicointerpolado(:, x_maxPosition))), picoapicointerpolado(:, x_maxPosition)/max(max(picoapicointerpolado)), 'o')
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes-50, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])

figure(90)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 60))), picoapico(:, 60)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = 0 mm')
xlabel('z (mm)')
ylabel('Normalizado')

figure(91)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 59))), picoapico(:, 59)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes - 5, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = -0,5 mm')
xlabel('z (mm)')
ylabel('Normalizado')

figure(92)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 58))), picoapico(:, 58)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes - 10, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = -1,0 mm')
xlabel('z (mm)')
ylabel('Normalizado')

figure(93)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 57))), picoapico(:, 57)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes - 15, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = -1,5 mm')
xlabel('z (mm)')
ylabel('Normalizado')

figure(94)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 56))), picoapico(:, 56)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes - 20, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = -2,0 mm')
xlabel('z (mm)')
ylabel('Normalizado')

figure(95)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 55))), picoapico(:, 55)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes - 25, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = -2,5 mm')
xlabel('z (mm)')
ylabel('Normalizado')

figure(96)
% Experimental signal - Pure
hold on; plot(0.5*(1:length(picoapico(:, 50))), picoapico(:, 50)/max(max(picoapico))); 
% Simulated signal
hold on; plot(z*1e3, abs(Ppp(rowsOfMaxes - 50, :))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
legend('Experimental signal', 'Simulated signal')
title ('Foco = 40 mm e x = -5,0 mm')
xlabel('z (mm)')
ylabel('Normalizado')


% Velocidade normal (input)
vn = resample(hdf5read('ref_pulse-40MHz.h5', 'ascan'), 2*f*sample*1e-6, 40);
plot(1e+6/(2*f*sample)*(1:length(vn)), vn)
grid on
grid minor
axis square
xlabel('Tempo $$({\mu s})$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('Amplitude Normalizada', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Sinal de excitação experimental', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
legend('vn' )
set(gca,'FontSize',20);