function [xscanned, yscanned, zscanned, Nscanned, xu] = scanner(xc, yc, zc, D, discretization) % Main function

% [xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization)
% xc --> x-coordinate of the center of the hydrophone
% yc --> y-coordinate of the center of the hydrophone
% zc --> z-coordinate for translation to the coordinate xy-plane
% D --> Diameter of the hydrophone
% discretization --> Value to divide the geometry into finite elements to prepare for analysis

R = D/2; % Radius of the hydrophone
xscan = linspace(xc - R, xc + R, discretization);
yscan = linspace(yc - R, yc + R, discretization);

[Xscan, Yscan] = meshgrid(xscan, yscan);  % Define um conjunto de pontos numa malha retangular.
Zscan = zc*ones(size(Yscan));
xscan = reshape(Xscan, [prod(size(Xscan))  1]); % Transforma em vetor as coordenadas x.
yscan = reshape(Yscan, [prod(size(Yscan))  1]); % Transforma em vetor as coordenadas y.
zscan = reshape(Zscan, [prod(size(Zscan))  1]); % Transforma em vetor as coordenadas z.

% Plota os pontos da área quadrada
% plot3(zscan, xscan, yscan, 'b.'); grid on; hold on;

f = (xscan - xc).^2 + (yscan - yc).^2 - R^2;
index = (f <= 0); % Pontos no interior com valor lógico '1'. Pontos fora = '0'.
xscanned = xscan(index); % Copia apenas as coordenadas dos pontos que estão dentro da área circular.
yscanned = yscan(index);
zscanned = zscan(index);
Nscanned = length(xscanned); % Número de pontos dentro da região circular de raio R.

% Identifica os pontos do eixo x no primeiro quadrante do hidrofone
xu = unique(xscanned);
yu = zeros(length(xu), 1);
zu = zc*ones(length(xu), 1);

% Plota os pontos da área circular
% plot3(zscanned, xscanned, yscanned, 'ro'); 

% Plota os pontos do eixo x no primeiro quadrante do hidrofone
% plot3(zu, xu, yu, 'b*'); hold off;