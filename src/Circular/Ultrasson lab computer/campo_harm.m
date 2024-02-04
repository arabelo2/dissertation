%function [xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization) % Main function

% [xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization)
% xc --> x-coordinate of the center of the hydrophone
% yc --> y-coordinate of the center of the hydrophone
% zc --> z-coordinate for translation to the coordinate xy-plane
% D --> Diameter of the hydrophone
% discretization --> Value to divide the geometry into finite elements to prepare for analysis

% clear all;
% close all;

D = 0.01905; % Diameter [m]

discretization = 100;

% Tx na origem
xc = 0;
yc = 0;
zc = 0; 

R = D/2; % Radius [m]
c1 = 1500; % [m/s]
rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
fs = 64e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
STEP = lambda / 10;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi/lambda; % wave number

xscan = linspace(xc - R, xc + R, discretization);
yscan = linspace(yc - R, yc + R, discretization);

[Xscan, Yscan] = meshgrid(xscan, yscan);  % Define um conjunto de pontos numa malha retangular.

Zscan = zc*ones(size(Yscan));
xscan = reshape(Xscan, [prod(size(Xscan))  1]); % Transforma em vetor as coordenadas x.
yscan = reshape(Yscan, [prod(size(Yscan))  1]); % Transforma em vetor as coordenadas y.
zscan = reshape(Zscan, [prod(size(Zscan))  1]); % Transforma em vetor as coordenadas z.

% Plota os pontos da área quadrada
plot3(zscan, xscan, yscan, 'b.'); grid on; hold on;

f = (xscan - xc).^2 + (yscan - yc).^2 - R^2;
index = (f <= 0); % Pontos no interior com valor lógico '1'. Pontos fora = '0'.
xscanned = xscan(index); % Copia apenas as coordenadas dos pontos que estão dentro da área circular.
yscanned = yscan(index);
zscanned = zscan(index);
Nscanned = length(xscanned); % Número de pontos dentro da região circular de raio R.

% Plota os pontos da área circular
plot3(zscanned, xscanned, yscanned, 'ro'); hold off;

%%%%%%%%%%%%

zmin = 0.002;  % z-axis
zmax = 0.102; % z-axis

xmin = -0.030; % x-axis
xmax = +0.030; % x-axis

zpoints = 667;
xpoints = 407;

z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);

r = cell(length(x), length(z));
pressure = cell(length(x), length(z));
Pp_c = zeros(length(x), length(z));

   
%sigma = sqrt((xscanned).^2 + (yscanned).^2 + (zscanned).^2);
%r = sqrt(x.^2 + 0.^2 + z.^2);

for zz = 1:length(z)
    for xx = 1:length(x)
        r{xx, zz} = sqrt((xscanned-x(xx)).^2 + (yscanned-0).^2 + (zscanned - z(zz)).^2);
        pressure{xx, zz} = sum(exp(-j*k*r{xx, zz})./r{xx, zz});

        % Peak amplitude
        Pp_c(xx, zz) = abs(pressure{xx, zz}); 
    end
end

Pp_c_analytics = Pp_c;

figure()
% plot(z*lambda/(R^2), Pp_c_analytics(61,:)/max(Pp_c_analytics(61,:)),'r')
plot(z*lambda/(R^2), Pp_c_analytics(floor(length(x)/2) + 1, :)/max(max(Pp_c_analytics)))

ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('$$z \frac{\lambda}{a^2}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressão ao longo do eixo acústico (eixo z)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor

hold on
Pp_opb = abs(exp(-j*k*z) - exp(-j*k*sqrt(z.^2 + R^2)));
plot(z*lambda/(R^2), Pp_opb/max(Pp_opb),'ko')

legend('Pp__c(61,:) (Integral de Rayleigh)','Pp__opb (Ondas Planas + Bordas')

figure
surface(Pp_c_analytics)
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor


pcolor(z, x, Pp_c_analytics/max(max(Pp_c_analytics)))
shading interp
colormap(jet)
colorbar('eastoutside')
axis padded
xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
zlabel('Pressão Normalizada', 'Color', 'k', 'interpreter', 'latex')
title('Pressao gerada através da Integral de Rayleigh', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect ('auto') % daspect([1 1 1])

% figure()
% pcolor(z, x, Pp_c/max(Pp_c))
% xlabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
% ylabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
% title(['Pressure field of a circular plane rigid baffled piston'], 'Color', 'k', 'interpreter', 'latex')
% az = 0; % az = -90/90 -- > Horizontal ; % az = 0/180 -- > Vertical;
% el = 90;
% view(az, el);
% shading interp
% colormap(jet)
% colorbar
% axis padded
% grid on
% grid minor
% set(gca,'Ydir','reverse');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
% daspect('auto') % daspect([1 1 1])
% 
figure()
mesh(z, x, Pp_c/max(max(Pp_c)))
xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
zlabel('Normalized pressure', 'Color', 'k', 'interpreter', 'latex')
title(['Pressure field of a circular plane rigid baffled piston through Rayleigh integral'], 'Color', 'k', 'interpreter', 'latex')
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect('auto') % daspect([1 1 1])

figure()
plot(z*lambda/(R^2), Pp_c_prg.Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg.Pp_c_prg)), 'k', z*lambda/(R^2), Pp_c_analytics(floor(length(x)/2) + 1, :)/max(max(Pp_c_analytics)), 'r')
ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('$$z \frac{\lambda}{a^2}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressão ao longo do eixo acústico (eixo z)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
legend('Pressão gerada através do programa','Pressão gerada através da Integral de Rayleigh')

figure ()
plot(x/R, Pp_c_prg.Pp_c_prg(:, floor(length(z)/2) + 58)/max(max(Pp_c_prg.Pp_c_prg)), 'k', x/R, Pp_c_analytics(:, floor(length(z)/2) + 58)/max(max(Pp_c_analytics)), 'r')
ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('Raio normalizado = $$\frac{x}{R}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Perfil de pressão acústica (eixo x) (z ~= NF) (NF = 0.06048375 m)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
legend('Pressão gerada através do programa','Pressão gerada através da Integral de Rayleigh')

figure ()
plot(x/R, Pp_c_prg.Pp_c_prg(:, 52)/max(max(Pp_c_prg.Pp_c_prg)), 'k', x/R, Pp_c_analytics(:, 52)/max(max(Pp_c_analytics)), 'r')
ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('Raio normalizado = $$\frac{x}{R}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Perfil de pressão acústica (eixo x) (z ~= R) (D/2 = R = 0.009525 m)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
legend('Pressão gerada através do programa','Pressão gerada através da Integral de Rayleigh')
