Ppp = load('fs40'); Ppp = Ppp.Ppp;
[rowsOfMaxes colsOfMaxes] = find(Ppp == max(max(Ppp)));
xMin = 2;
xMax = 102;
yMin = -30;
yMax = 30;

figure(50)
hold on; plot(z(length(find(z < xMin*1e-3))+1:end)*1e3, abs(Ppp(rowsOfMaxes, length(find(z < xMin*1e-3))+1:end))./max(max(Ppp)))
grid on
grid minor
ylim([0,1])
% legend('fs = 640MHz', 'fs = 320MHz', 'fs = 160MHz', 'fs = 100MHz', 'fs =  80MHz', 'fs =  40MHz')
% legend('fs temporal = 160MHz e fs espacial =  13.33 ciclos por mm', 'fs = 100MHz', 'fs =  40MHz')
xlabel('z (mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel(' Amplitude Normalizado', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title ('Pressão acústica ao longo do eixo')
set(gca,'FontSize',20);