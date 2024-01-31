%function [xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization) % Main function

% [xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization)
% xc --> x-coordinate of the center of the hydrophone
% yc --> y-coordinate of the center of the hydrophone
% zc --> z-coordinate for translation to the coordinate xy-plane
% D --> Diameter of the hydrophone
% discretization --> Value to divide the geometry into finite elements to prepare for analysis

clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 200;
xc = 0; yc = 0; zc = 0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
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
%Nscanned = length(xscanned); % Número de pontos dentro da região circular de raio R.

% Plota os pontos da área circular
plot3(zscanned, xscanned, yscanned, 'ro'); hold off;

%%%%%%%%%%%%

zmin = 0.002;  % z-axis
zmax = +0.102; % z-axis

xmin = -0.030; % x-axis
xmax = +0.030; % x-axis

zpoints = 201; 
xpoints = 121;

z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);

for zz = 1:length(z)
    for xx = 1:length(x)    

        r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
        pressao = sum(exp(-j*k*r)./r);
        campo(xx, zz) = pressao;
               
        % Peak amplitude
        Pp_c(xx, zz) = abs(pressao); 
    end
end

Pp_opb = abs(exp(-j*k*z) - exp(-j*k*sqrt(z.^2 + R^2)));

figure
surface(Pp_c)
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor

figure()
pcolor(z, x, Pp_c/max(Pp_c))
xlabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
title(['Pressure field of a circular plane rigid baffled piston'], 'Color', 'k', 'interpreter', 'latex')
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
% 
figure()
mesh(z, x, Pp_c/max(Pp_c))
xlabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
% zlabel('Max(|P(r,t)|)')
zlabel('Normalized pressure', 'Color', 'k', 'interpreter', 'latex')
title(['Pressure field of a circular plane rigid baffled piston '], 'Color', 'k', 'interpreter', 'latex')
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect('auto') % daspect([1 1 1])
