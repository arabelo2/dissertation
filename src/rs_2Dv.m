function p = rs_2Dv(b, f, c, e, x, z, varargin)
% Compute wave number
kb = 2*pi.*b.*f/c;
% If a number of segments is specified, use:
if nargin == 7
    N = varargin{1};
else
    N = round(20*f*b/c);
    if N < 1
        N = 1;
    end
end
% Use normalized positions in the fluid.
xb = x./b;
zb = z./b;
eb = e./b;
% Compute normalized centroid locations for the segments.
xc = zeros(1, N);
for jj = 1:N
    xc(jj) = -1 + 2*(jj - 0.5)/N;
end
% Calculate normalized pressure as a sum over all the segments as an
% approximation of the Rayleigh-Somerfeld type of integral.
p = 0;
for kk = 1:N
    rb = sqrt((xb - xc(kk) - eb).^2 + zb.^2);
    p  = p + besselh(0, 1, kb.*rb);
end
p = p.*(kb./N); % Include external factor

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.