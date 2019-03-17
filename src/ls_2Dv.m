function p = ls_2Dv(b, f, c, e, x, z, varargin)
% Compute wave number
kb = 2*pi.*b.*f/c;
% If number of segments is specified, use:
if nargin == 7
    N = varargin{1};
else
    N = round(20*f*b/c);
    if N < 1
        N = 1;
    end
end
% Use normalized positions in the fluid
xb = x./b;
zb = z./b;
eb = e./b;
% Compute normalized centroid locations for the segments
xc = zeros(1,N);
for jj = 1:N
    xc(jj) = -1 + 2*(jj - 0.5)/N;
end
% Calculate normalized pressure as a sum over all the segments as an
% approximation of the Rayleigh-Sommerfeld tye of integral
p = 0;
for kk = 1:N
    ang = atan((xb - xc(kk) - eb)./zb);
    ang = ang + eps.*(ang == 0);
    dir = sin(kb.*sin(ang)/N)./(kb.*sin(ang)/N);
    rb = sqrt((xb - xc(kk) - eb).^2 + zb.^2);
    ph = exp(i*kb.*rb);
    p = p + dir.*exp(i*kb.*rb)./sqrt(rb);
end
p = p.*(sqrt(2*kb./(i*pi)))/N; % include external factor

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.