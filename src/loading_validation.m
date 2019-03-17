tic
f0 = 5e6;
c = 1480;
c1 = c;
rho1  = 1000;
c2 = 1480;
rho2  = 1000;
b = 6e-3/2;
e = 0*b;
lambda = c/f0;
Dt0 = 50.8e-3;
angt = 10.217;

d1 = 1000;   % Density, medium one (arbitrary units)
cp1 = 1480; % Compressional wave speed, medium one (m/s)
d2 = 7900;   % Density, medium two (arbitrary units)
cp2 = 5900; % Compressional wave speed, medium two (m/s)
cs2 = 3200; % Shear wave speed, medium two (m/s)
type = 'p'; % Wave type, medium two
mat = {d1, cp1, d2, cp2, cs2, type};    % Form material vector

% Field points (x, y, z) to evaluate
N = 10;
xmin = 0e-3;
xmax = 25e-3;
xnpoints = N*ceil(abs(xmax - xmin)/lambda);
xs = linspace(xmin, xmax, xnpoints);
zmin = 1e-3;
zmax = 25e-3;
znpoints = N*ceil(abs(zmax - zmin)/lambda);
zs = linspace(zmin, zmax, znpoints);
y = 0;
[x, z] = meshgrid(xs, zs);

if 2*b > lambda/10
    Nopt = ceil(20*f0*b/c);
else
    Nopt = 1;
end

figure(1)
p_ls_2Dint = ls_2Dint(b, f0, mat, e, angt, Dt0, x, z, Nopt);
imagesc(1000*xs, 1000*zs, abs(p_ls_2Dint))
shading interp
colormap(gray)
colorbar
axis vis3d
xlabel('x, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('z, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressao normalizada (ls_2Dint)',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);

% ----------------------------------------------------------------------

% figure(1)
% p_rs = rs_2Dv(b, f0, c, e, x, z, Nopt);
% imagesc(1000*xs, 1000*zs, abs(p_rs))
% shading interp
% colormap(gray)
% colorbar
% axis vis3d
% xlabel('x, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('z, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressao normalizada (rs)',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% ----------------------------------------------------------------------

% figure(2)
% p_ls = ls_2Dv(b, f0, c, e, x, z, Nopt);
% imagesc(1000*xs, 1000*zs, abs(p_ls))
% shading interp
% colormap(gray)
% colorbar
% axis vis3d
% xlabel('x, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('z, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressao normalizada (ls)',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% ----------------------------------------------------------------------

% xs = 0;
% p = rs_2Dv(b, f0, c, e, xs, zs, Nopt);
% plot(1000*zs, abs(p), 'b')
% hold on
% plot(1000*zs, abs(p), 'ko')
% hold off
% axis normal
% xlabel('z, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('$$|\frac{p}{\rho c v}|$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressao normalizada sobre o eixo',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% ----------------------------------------------------------------------

% zs = 100e-3;
% p = rs_2Dv(b, f0, c, e, xs, zs, Nopt);
% plot(rad2deg(asin(xs/zs)), abs(p)/max(abs(p)))
% hold on
% plot(rad2deg(asin(xs/zs)), abs(p)/max(abs(p)), 'ko')
% hold off
% axis normal
% ylim([0 1])
% xlabel('$$\theta, degrees$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('$$|\frac{p(r_0, \theta)}{p(r_0, 0)}|$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressao normalizada sobre o eixo',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

toc