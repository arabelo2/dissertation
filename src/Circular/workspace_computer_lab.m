loadingDataOfRectangularArrayPistonWScanner
pwd
dir
cd src
dir
loadingDataOfRectangularArrayPistonWScanner
dir C:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array
cd 'C:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array'
loadingDataOfRectangularArrayPistonWScanner
%-- 29/01/2024 13:07 --%
run_c
.1/STEP
.061/STEP
run_c
run_c_wscanner
size(h_c)
discretization = 40;
D = 19e-3;
xc = 0; yc = 0; % Tx na origem
R = D/2; % Radius of transducer
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
discretization = 40;
D = 19e-3;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius of transducer
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
clear all, clc
campo_harm
clear all, clc
D = 0.019; % Diameter [m]
discretization = 40;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xmin = 0.002;  % z-axis
xmax = +0.102; % z-axis
ymin = -0.030; % x-axis
ymax = +0.030; % x-axis
xpoints = 201;
ypoints = 120;
x = linspace(xmin, xmax, xpoints);
y = linspace(ymin, ymax, ypoints);
for xx = 1:length(x)
for yy = 1:length(y)
r = sqrt((xscan-x(xx)).^2 + (yscan-y(yy)).^2);
pressao(yy, xx) = sum(exp(-j*k*r)./r);
% Peak amplitude
Pp_c(yy, xx) = abs(pressao);
end
end
clear
clear all
D = 0.019; % Diameter [m]
discretization = 40;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xmin = 0.002;  % z-axis
xmax = +0.102; % z-axis
ymin = -0.030; % x-axis
ymax = +0.030; % x-axis
xpoints = 201;
ypoints = 120;
x = linspace(xmin, xmax, xpoints);
y = linspace(ymin, ymax, ypoints);
xx = 1
yy = 1
r = sqrt((xscan-x(xx)).^2 + (yscan-y(yy)).^2);
pressao(yy, xx) = sum(exp(-j*k*r)./r);
clear all
campo_harm
clear all
campo_harm
clear all
campo_harm
clear all
D = 0.019; % Diameter [m]
discretization = 40;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xmin = 0.002;  % z-axis
xmax = +0.102; % z-axis
ymin = -0.030; % x-axis
ymax = +0.030; % x-axis
xpoints = 201;
ypoints = 120;
x = linspace(xmin, xmax, xpoints);
y = linspace(ymin, ymax, ypoints);
for xx = 1:length(x)
for yy = 1:length(y)
r = sqrt((xscan-x(xx)).^2 + (yscan-y(yy)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(yy, xx) = pressao;
% Peak amplitude
Pp_c(yy, xx) = abs(pressao);
end
end
figure
surface(Pp_c)
j
% Plota os pontos da área quadrada
plot3(zscan, xscan, yscan, 'b.'); grid on; hold on;
% Plota os pontos da área circular
plot3(zscanned, xscanned, yscanned, 'ro'); hold off;
clear all
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
zmin = 0.002;  % z-axis
zmax = +0.102; % z-axis
xmin = -0.030; % x-axis
xmax = +0.030; % x-axis
zpoints = 201;
xpoints = 120;
z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);
for zz = 1:length(z)
for xx = 1:length(x)
r = sqrt((xscan-x(xx)).^2 + (yscan-0)).^2 + (zscan - z(zz).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
% Peak amplitude
Pp_c(xx, zz) = abs(pressao);
end
end
figure
surface(Pp_c)
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
size(Pp_c)
size(x)
figure()
pcolor(x, z, Pp_c/max(Pp_c))
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
surface(Pp_c)
doc surface
[X,Z] = meshgrid(x, z);
surface(X,Z, Pp_c)
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
zmin = 0.002;  % z-axis
zmax = +0.102; % z-axis
xmin = -0.030; % x-axis
xmax = +0.030; % x-axis
zpoints = 201;
xpoints = 120;
z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);
x(1)
zz = 1
xx = 1
xscan-x(xx)
yscan-0
r = sqrt((xscan-x(xx)).^2 + (yscan-0)).^2 + (zscan - z(zz)).^2);
clear
close all
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xpoints = 120;
z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);
zz = 1
xx = 1
r = sqrt((xscan-x(xx)).^2 + (yscan-0)).^2 + (zscan - z(zz)).^2);
(zscan - z(zz))
(zscan - z(zz)).^2
size((zscan - z(zz)))
size((xscan-x(xx)))
size((yscan-0))
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2;
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
Pp_c(xx, zz) = abs(pressao);
campo_harm
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xpoints = 120;
z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);
for zz = 1:length(z)
for xx = 1:length(x)
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2;
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
% Peak amplitude
Pp_c(xx, zz) = abs(pressao);
end
end
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xpoints = 120;
z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);
zz = 1
xx = 1
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2;
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
Pp_c(xx, zz) = abs(pressao)
pressao
pressao = sum(exp(-j*k*r)./r);
pressao
Pp_c(xx, zz) = abs(pressao)
size(Pp_c)
zz =2
xx =2
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2;
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
Pp_c(xx, zz) = abs(pressao);
size(Pp_c)
Pp_c
zz = 1
xx = 1
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
Pp_c{xx, zz} = abs(pressao);
zz = 1
xx = 1
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
Pp_c(xx, zz) = abs(pressao)
clear all
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xpoints = 120;
z = linspace(zmin, zmax, zpoints);
x = linspace(xmin, xmax, xpoints);
zz = 1
xx = 1
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
% Peak amplitude
Pp_c(xx, zz) = abs(pressao);
Pp_c
xx =2
r = sqrt((xscan-x(xx)).^2 + (yscan-0).^2 + (zscan - z(zz)).^2);
pressao = sum(exp(-j*k*r)./r);
campo(xx, zz) = pressao;
Pp_c(xx, zz) = abs(pressao);
Pp_c
clear all;
close all;
D = 0.019; % Diameter [m]
discretization = 100;
xc = 0; yc = 0; zc=0; % Tx na origem
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
%rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
%fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
%STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
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
xpoints = 120;
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
surface(Pp_c)
campo_harm
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
figure
surface(Pp_c)
shading interp
colormap(jet)
colorbar
campo_harm
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
zscan(1:3)
xscan(1:10)
yscan(1:10)'
campo_harm
x(60)
x(61)
campo_harm
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
figure
surface(Pp_c)
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
x(61)
figure
plot(z,Pp_c(61,:))
k
k*R
p300=Pp_c(61,:);
campo_harm
figure
plot(p300)
hold on
p400=Pp_c(61,:);
plot(p400,'r')
p_ondapb = exp(-j*k*z) - exp(-j*k*sqrt(z.^2 + R^2));
p_ondapb
campo_harm
plot(z,Pp_opb)
run_o
plot(z,Pp_c(61, :))
plot(z,Pp_c(61, :)//max(max(Pp_c)))
plot(z,Pp_c(61, :)/max(max(Pp_c)))
Pp_opb = abs(exp(-j*k*z) - exp(-j*k*sqrt(z.^2 + R^2)));
hold on
plot(z,Pp_opb/max(Pp_opb),'.r')
dt=t_temp(2)-t_temp(1)
c1*dt
c1*dt*e3
c1*dt*1e3
plot(z,Pp_c(61, :)/max(max(Pp_c)), 'ko')