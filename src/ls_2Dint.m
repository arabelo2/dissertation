function p = ls_2Dint(b, f, mat, e, angt, Dt0, x, z, varargin)
% Extract material parameters
d1 = mat{1}; % Density in the first medium kg/m^3
c1 = mat{2}; % Wave speed in the first medium m/s
d2 = mat{3}; % Density in the second medium kg/m^3
c2 = mat{4}; % Wave speed in the second medium m/s
% Compute wave numbers
k1b = 2*pi.*b.*f/c1;
k2b = 2*pi.*b.*f/c2;
% If number of segments is specified, use:
if nargin == 9
    N = varargin{1};
else
    N = round(20*f*b/c1);
    if N < 1
        N = 1;
    end
end
% Compute centroid locations for the segments
xc = zeros (1, N);
for jj = 1:N
    xc(jj) = b*(-1 + 2*(jj - 0.5)/N);
end
% Calculate normalized pressure as a sum over all the segments
p = 0;
for nn = 1:N
    % Find the distance, xi, where the ray from the center of a segment
    % to point(x, z) intersects the interface
    xi = pts_2Dintf(e, xc(nn), angt, Dt0, c1, c2, x, z);
    % Compute the distances and angles needed in the model
    Dtn = Dt0 + (e + xc(nn)).*sin(angt*pi/180);
    Dxn = x - (e + xc(nn)).*cos(angt*pi/180);
    r1 = sqrt(xi.^2 + Dtn.^2)./b;
    r2 = sqrt((Dxn - xi).^2 + z.^2)./b;
    ang1 = asin(xi./(b*r1));
    ang2 = asin((Dxn - xi)./(b*r2));
    ang = angt*pi/180 - ang1;
    ang = ang + eps.*(ang == 0);
    % Form up the segment directivity
    dir = sin(k1b.*sin(ang)/N)./(k1b.*sin(ang)/N);
    % Compute plane wave transmission coefficient (based on pressure ratio)
    Tp = 2.*d2.*c2.*cos(ang1)./(d1.*c1.*cos(ang2) + d2.*c2.*cos(ang1));
    % Compute phase term and denominator
    ph = exp(1i*k1b.*r1 + 1i*k2b.*r2);
    den = r1 + (c2/c1).*r2.*((cos(ang1)).^2)./(cos(ang2)).^2;
    % Put terms together for pressure due to each segment
    p = p + Tp.*dir.*ph./sqrt(den);
end
p = p.*(sqrt(2*k1b./(1i*pi)))/N; % Include external factor
    
% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.