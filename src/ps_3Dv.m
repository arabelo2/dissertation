function p = ps_3Dv (lx, ly, f, c, ex, ey, x, y, z, varargin)
% Calculate wave number
k = 2*pi*f/c;
% If number of x-segments is specified then use
if nargin > 9
    P = varargin{1};
    % else choose number of terms so each segment length is at most a wave
    % length.
else
    P = ceil(f*lx/c);
    if P < 1
        P = 1;
    end
end
% If number of y-segments is specified then use
if nargin > 10
    Q = varargin{2};
    % else choose number of terms so that each segment is a wave length or
    % less.
else
    Q = ceil(f*ly/c);
    if Q < 1
        Q = 1;
    end
end
% Compute centroid locations of segments in x- and y- directions.
xc = zeros(1, P);
yc = zeros(1, Q);
for pp = 1:P
    xc(pp) = -lx/2 + lx/P*(pp - 0.5);
end
for qq = 1:Q
    yc(qq) = -ly/2 + ly/Q*(qq - 0.5);
end
% Calculate normalized pressure as a sum over all the segments as an
% approximation of the Rayleigh-Sommerfeld integral.
p = 0;
for pp = 1:P
    for qq = 1:Q
        rpq = sqrt((x - xc(pp) - ex).^2 + (y - yc(qq) - ey).^2 + z.^2);
        ux = (x - xc(pp) - ex)./rpq;
        uy = (y - yc(qq) - ey)./rpq;
        ux = ux + eps*(ux == 0);
        uy = uy + eps*(uy == 0);
        dirx = sin(k.*ux.*lx/(2*P))./(k.*ux.*lx/(2*P));
        diry = sin(k.*uy.*ly/(2*Q))./(k.*uy.*ly/(2*Q));
        p = p + dirx.*diry.*exp(1i*k.*rpq)./rpq;
    end
end
p = p.*(-1i*k*lx/P*ly/Q)/(2*pi); % Include external factor

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.