function p = rs_2Dv(b, f, c, e, x, z, varargin)
% p = rs_2Dv(b, f, c, e, x, z, Nopt) computes the normalized pressure, p,
% at a location (x, z) (in mm) in a fluid for a 1-D element of length 2b
% (in mm) along the x-axis at a frequency, f, (in MHz). And for a wave
% speed, c, (in m/sec) of the fluid. This function can be used to describe
% an element in an array by specifying a non-zero value, e, (in mm) which
% in the offset of the center of the element along the x-axis. The assumed
% harmonic time dependency is exp(-2i*pi*f*t) and the 2-D version of the
% Rayleigh-Sommerfeld integral for a piston source is used as the model.
% Nopt gives the number of segments to use. If Nopt is not given as an
% input argument the function uses 10 segments per wavelength, based on the
% input frequency, f, which must be a scalar when Nopt is not given.

% Compute wave number
kb = 2000*pi.*b.*f/c;

% If a number of segments is specified, use
if nargin == 7
    N = varargin{1};
else
    % Else choose number of segments so that the size of each segment is
    % one-tenth a wavelength.
    N = round(20000*f*b/c);
    if N <= 1
        N = 1;
    end
end

% Use normalized positions in the fluid
xb = x./b;
zb = z./b;
eb = e./b;

% Compute normalized centroid locations for the segments
xc = zeros(1, N);
for jj = 1:N
    xc(jj) = -1 + 2*(jj - 0.5)/N;
end

% Calculate normalized pressure as a sum over all the segments as an
% approximation of the Rayleigh-Somerfeld type of integral

p = 0;
%P = zeros(1, N);
for kk = 1:N
    rb = sqrt((xb - xc(kk) - eb).^2 + zb.^2);
    p  = p + besselh(0, 1, kb.*rb);
end

p = p.*(kb./N); % Include external factor

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.