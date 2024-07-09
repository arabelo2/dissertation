T = length(V{1, 1});
signals_matrix = zeros(length(x), length(z), T);
for zz = 1:length(z)
    for xx = 1:length(x)
        signals_matrix(xx, zz, :) = V{xx, zz}; % Each cell contains a vector of length T
    end
end

xxmin = -0.030; % x-axis
xxmax = 0.030; % x-axis
zzmin = 0; % z-axis
zzmax = 0.120; % z-axis
xxpoints = 121; xxx = floor(xxpoints/2) + 1;
zzpoints = 121; zzz = floor(zzpoints/2) + 1;
xq = linspace(xxmin, xxmax, xxpoints);
zq = linspace(zzmin, zzmax, zzpoints); zq = zq + 4.8375e-04;
[zq_grid, xq_grid] = meshgrid(zq, xq);
interpolated_signals = zeros(length(xq), length(zq), T);
[z_grid, x_grid] = meshgrid(z, x);	

rr = 1;
	
for t = 1:T
    interpolated_signals(:, :, t) = interp2(z_grid, x_grid, signals_matrix(:, :, t), zq_grid, xq_grid, 'linear');
    % interpolated_signals(:, :, t) = interp(signals_matrix(:, :, t), rr);
end

figure(1)
mesh(zq_grid, xq_grid, max(interpolated_signals, [], 3)/max(interpolated_signals, [], 'all'))

figure(2)
pcolor(zq_grid, xq_grid, max(interpolated_signals, [], 3)/max(interpolated_signals, [], 'all'))
xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
title(['Pressão Normalizada (Método Analítico) com', num2str(nc),' ciclo(s)'], 'Color', 'k', 'interpreter', 'latex')
az = 0; % az = -90/90 -- > Horizontal ; % az = 0/180 -- > Vertical;
el = 90;
view(az, el);
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect('auto') % daspect([1 1 1])

figure(10)
hold on
plot(zq*lambda/(R^2), max(interpolated_signals(xxx,:,:), [], 3)/max(interpolated_signals, [], 'all'))

figure(4)
hold on
plot(xq, max(interpolated_signals(:,zzz,:), [], 3)/max(interpolated_signals, [], 'all'))
