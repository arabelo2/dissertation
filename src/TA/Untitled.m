figure(1)
subplot(1,2,1)
imagesc(x_axis(x_idx), z_axis(z_idx), picoapico(z_idx, x_idx)/max(max(picoapico)))
colormap(jet)
colorbar
axis normal
ylabel('z(mm)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('x(mm)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title({'Amplitude normalizada pico-a-pico', 'Campo de pressão experimental'}, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);
daspect([1 1 1])

figure(1)
subplot(1,2,2)
imagesc(x*1e3, z*1e3, abs(Ppp')/max(max(Ppp)))
shading interp
colormap(jet)
colorbar
axis normal
ylabel('z(mm)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('x(mm)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title({'Amplitude normalizada pico-a-pico', 'Campo de pressão simulado'}, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca,'FontSize',20);
daspect([1 1 1])

figure()
hold on
% pe_l = plot(z_axis, picoapico(:,xe)/max(max(picoapico)), 'k');
pe_o = plot(z_axis, picoapico(:,xe)/max(max(picoapico)), 'b--o');
pa_l = plot(linspace(z_axis(1),z_axis(end), size(Ppp', 1)), Ppp(xa, :)'/max(max(Ppp)), 'r');
xlabel('z(mm)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('Pressão normalizada', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title({'Amplitude da pressão acústica ao longo do eixo'}, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);
hold off
legend([pe_o pa_l], 'Experimental', 'Simulada')

figure()
hold on
%pe_l = plot(x_axis, picoapico(ze,:)/max(max(picoapico)), 'k--o');
pe_o = plot(x_axis, picoapico(ze,:)/max(max(picoapico)), 'b--o');
pa_l = plot(linspace(x_axis(1),x_axis(end), size(Ppp', 2)), Ppp(:, za)'/max(max(Ppp)), 'r');
xlabel('x(mm)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('Pressão normalizada', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title({'Perfil da pressão acústica' , 'na posição do foco'}, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);
hold off
legend([pe_o pa_l], 'Experimental', 'Simulado')


figure()
hold on
%pe_l = plot(x_axis, picoapico(ze,:)/max(max(picoapico)), 'k--o');
pe_o = plot(rad2deg(tan(x_axis/x_axis(end))), mag2db(picoapico(ze,:)/max(max(picoapico))), 'b--o');
pa_l = plot(rad2deg(tan(linspace(x_axis(1), x_axis(end), size(Ppp', 2))/x_axis(end))), mag2db(Ppp(:, za)'/max(max(Ppp))), 'r');
xlabel('${\theta ^{\circ}}$', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('Magnitude (dB)', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title({'Perfil da pressão acústica' , 'na posição do foco'}, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);
hold off
legend([pe_o pa_l], 'Experimental', 'Simulado')


p1 = plot(mean(B));
set(p1, 'LineWidth', 2, 'MarkerSize', 10, 'Marker', 'x')
hold on
e1 = errorbar(mean(B), sem);
set(e1, 'LineStyle', 'none');
eline = get(e1, 'Children');
set(eline,  'LineWidth', 2, 'Marker', 'x', 'MarkerSize', 4)

% With mag2db
figure()
hold on
ordem1 = 43
phase1 = rad2deg(tan(x_axis/x_axis(end)))
mag1 = mag2db(picoapico(ze,:)/max(max(picoapico)))
plot(phase1, mag1, 'r--o')
p1 = polyfit(phase1, mag1, ordem1);
plot(phase1, polyval(p1, phase1), 'k')
[p1 s1] = polyfit(phase1, mag1, ordem1);
[y_fit1 delta1] = polyval(p1, phase1, s1);
plot(phase1, y_fit1 + 2 * delta1, 'k:', phase1, y_fit1 - 2 * delta1, 'k:')
phase2 = rad2deg(tan(linspace(x_axis(1), x_axis(end), size(Ppp', 2))/x_axis(end)))
mag2 =  mag2db(Ppp(:, za)'/max(max(Ppp)))
plot(phase2, mag2, 'b--.')
grid on
grid minor
hold off


figure()
hold on
ordem2 = 50;
phase2 = rad2deg(tan(linspace(x_axis(1), x_axis(end), size(Ppp', 2))/x_axis(end)))
mag2 =  mag2db(Ppp(:, za)'/max(max(Ppp)))

plot(phase2, mag2, '.')
p2 = polyfit(phase2, mag2, ordem2);
plot(phase2, polyval(p2, phase2))
[p2 s2] = polyfit(phase2, mag2, ordem2);
[y_fit2 delta2] = polyval(p2, phase2, s2);
plot(phase2, y_fit2 + 2 * delta2, ':', phase2, y_fit2 - 2 * delta2, ':')

phase1 = rad2deg(tan(x_axis/x_axis(end)))
mag1 = mag2db(picoapico(ze,:)/max(max(picoapico)))
plot(phase1, mag1, '*')
hold off

% Without mag2db
figure()
hold on
phase1 = rad2deg(tan(x_axis/x_axis(end)))
mag1 = (picoapico(ze,:)/max(max(picoapico)))
plot(phase1, mag1, '.')
p1 = polyfit(phase1, mag1, 5);
plot(phase1, polyval(p1, phase1))
[p1 s1] = polyfit(phase1, mag1, 5);
[y_fit1 delta1] = polyval(p1, phase1, s1);
plot(phase1, y_fit1 + 2 * delta1, ':', phase1, y_fit1 - 2 * delta1, ':')
phase2 = rad2deg(tan(linspace(x_axis(1), x_axis(end), size(Ppp', 2))/x_axis(end)))
mag2 =  (Ppp(:, za)'/max(max(Ppp)))
plot(phase2, mag2, '*')
hold off

figure()
hold on

phase2 = rad2deg(tan(linspace(x_axis(1), x_axis(end), size(Ppp', 2))/x_axis(end)))
mag2 =  (Ppp(:, za)'/max(max(Ppp)))

plot(phase2, mag2, '.')
p2 = polyfit(phase2, mag2, 5);
plot(phase2, polyval(p2, phase2))
[p2 s2] = polyfit(phase2, mag2, 5);
[y_fit2 delta2] = polyval(p2, phase2, s2);
plot(phase2, y_fit2 + 2 * delta2, ':', phase2, y_fit2 - 2 * delta2, ':')

phase1 = rad2deg(tan(x_axis/x_axis(end)))
mag1 = (picoapico(ze,:)/max(max(picoapico)))
plot(phase1, mag1, '*')
hold off