function [xscanned, yscanned, zscanned, ru] = scanner(xc, yc, zc, D, discretization) % Main function

% [xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization)
% xc --> x-coordinate of the center of the hydrophone
% yc --> y-coordinate of the center of the hydrophone
% zc --> z-coordinate for translation to the coordinate xy-plane
% D --> Diameter of the hydrophone
% discretization --> Value to divide the geometry into finite elements to prepare for analysis

R = D/2; % Radius of the hydrophone

if discretization > 2
    xscan = linspace(xc - R, xc + R, discretization);
    yscan = linspace(yc - R, yc + R, discretization);
elseif (discretization <= 2) && (discretization > 0)
    xscan = xc;
    yscan = yc;
else
    disp("Error message: discretization must be greater than 0.")
    return
end

[Xscan, Yscan] = meshgrid(xscan, yscan);  % Define um conjunto de pontos numa malha retangular.
Zscan = zc*ones(size(Yscan));
xscan = reshape(Xscan, [prod(size(Xscan))  1]); % Transforma em vetor as coordenadas x.
yscan = reshape(Yscan, [prod(size(Yscan))  1]); % Transforma em vetor as coordenadas y.
zscan = reshape(Zscan, [prod(size(Zscan))  1]); % Transforma em vetor as coordenadas z.

% Plota os pontos da área quadrada
%plot3(zscan*1e3, xscan*1e3, yscan*1e3, 'k.'); grid on; hold on;

f = (xscan - xc).^2 + (yscan - yc).^2 - R^2;
index = (f <= 0); % Pontos no interior com valor lógico '1'. Pontos fora = '0'.
xscanned = xscan(index); % Copia apenas as coordenadas dos pontos que estão dentro da área circular.
yscanned = yscan(index);
zscanned = zscan(index);
% Nscanned = length(xscanned); % Número de pontos dentro da região circular de raio R.

% Identifica os pontos do eixo x do hidrofone
% xu = unique(xscanned);
% yu = zeros(length(xu), 1);
% zu = zc*ones(length(xu), 1);

% Calcula as distâncias de cada ponto no hidrofone em
ru = sqrt(xscanned.^2 + yscanned.^2);

% Plota os pontos da área circular
%plot3(zscanned*1e3, xscanned*1e3, yscanned*1e3, 'ro'); 

% Plota os pontos do eixo x no primeiro quadrante do hidrofone
%plot3(zscanned*1e3, xscanned*1e3, yscanned*1e3, 'b*');
% plot3(zu, xu, yu, 'k*'); 
%hold off;