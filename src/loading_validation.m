f0 = 5e6;
c = 1500;
b = .03e-3/2;
e = .20*b;

lambda = c/f0;

% Field points (x, y, z) to evaluate
N = 10;
xmin = -10e-3;
xmax = 10e-3;
xnpoints = N*ceil(abs(xmax - xmin)/lambda);
xs = linspace(xmin, xmax, xnpoints);

zmin = 0e-3;
zmax = 20e-3;
znpoints = N*ceil(abs(zmax - zmin)/lambda);
zs = linspace(zmin, zmax, znpoints);

y = 0;
[x, z] = meshgrid(xs, zs);

if 2*b > lambda/10
    Nopt = ceil(20*f0*b/c);
else
    Nopt = 1;
end

p = rs_2Dv(b, f0, c, e, x, z, Nopt);

imagesc(1000*xs, 1000*zs, abs(p))
shading interp
colormap(jet)
colorbar
axis vis3d
xlabel('x, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('z, mm', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressao normalizada',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);