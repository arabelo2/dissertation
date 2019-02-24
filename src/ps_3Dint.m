function [vx, vy, vz] = ps_3Dint(lx, ly, f, mat, ex, ey, angt, Dt0, x, y, z, varargin)
% Extract material densities, wave speeds, and the type of wave in the
% second medium from mat vector.
d1 = mat{1};
cp1 = mat{2};
d2 = mat{3};
cp2 = mat{4};
cs2 = mat{5};
type = mat{6};
% Wave speed in the first medium (a fluid) is for compressional waves
c1 = cp1;
% Decide which wave speed to use in second medium for specified wave type
if strcmp (type, 'p')
    c2 = cp2;
elseif strcmp (type, 's')
    c2 = cs2;
else
    error('Type must be ''p'' or ''s''')
end
% Compute wave numbers for waves in first and second medium
k1 = 2*pi*f/c1;
k2 = 2*pi*f/c2;
% If number of x-segments is specified then use
if nargin > 11
    R = varargin{1};
% Else choose number of terms so each segment is a wave length or less
else
    R = ceil(f*lx/c1);
    if R < 1
        R = 1;
    end
end
% If number of y-segments is specified then use
if nargin > 12
    Q = varargin{2};
% Else choose number of terms so that each segment is a wave length or less
else
    Q = ceil(f*ly/c1);
    if Q < 1
        Q = 1;
    end
end
% Compute centroid locations of segments in x' - and y' -directions
% relative to the element centroid
xc = zeros(1, R);
yc = zeros(1, Q);
for rr = 1:R
    xc(rr) = -lx/2 + lx/R*(rr - 0.5);
end
for qq = 1:Q
    yc(qq) = -ly/2 + ly/Q*(qq - 0.5);
end
% Calcule normalized velocity components as a sum over all the segments as
% an approximation of the Rayleight-Sommerfeld integral.
vx = 0;
vy = 0;
vz = 0;
for rr = 1:R
    for qq = 1:Q
        % Calculate distance xi along the interface for a ray from a
        % segment to the specified point in the second medium.
        Db = sqrt((x - (ex + xc(rr)).*cosd(angt)).^2 + (y - (ey + yc(qq))).^2);
        Ds = Dt0 + (ex + xc(rr)).*sind(angt);
        xi = pts_3Dint(ex, ey, xc(rr), yc(qq), angt, Dt0, c1, c2, x, y, z);
        % Calculate incident and refracted angles along the ray, including
        % the special case when ray is at normal incidence
        if Db == 0
            ang1 = 0;
        else
            ang1 = atand(xi./Ds);
        end
        if ang1 == 0
            ang2 = 0;
        else
            ang2 = atand((Db - xi)./z);
        end
        % Calculate ray path lengths in each medium
        r1 = sqrt(Ds.^2 + xi.^2);
        r2 = sqrt((Db - xi).^2 + z.^2);
        % Calculate segment sizes in x'- and y'- directions
        dx = lx/R;
        dy = ly/Q;
        % Calculate (x', y') components of unit vector along the ray in the
        % first medium
        if Db == 0
            uxt = -sind(angt);
            uyt = 0;
        else
            uxt = xi.*(x - (ex +xc(rr)).*cosd(angt)).*cosd(angt)./(Db.*r1) - Ds.*sind(angt)./r1;
            uyt = xi.*(y - (ey + yc(qq)))./(Db.*r1);            
        end
        % Calculate polarization components for P- and S-waves in the
        % second medium, including special case of normal incidence.
        if Db == 0
            dpx = 0;
            dpy = 0;
            dpz = 1;
            dsx = 1;
            dsy = 0;
            dsz = 0;
        else
            dpx = (1 - xi./Db).*(x - (ex + xc(rr)).*cosd(angt))./r2;
            dpy = (1 - xi./Db).*(y - (ey + yc(qq)))./r2;
            dpz = z./r2;
            dsx = sqrt(dpy.^2 + dpz.^2);
            dsy = -dpx.*dpy./dsx;
            dsz = -dpx.*dpz./dsx;
        end
        % Choose polarization components to use based on wave type in the
        % second medium
        if strcmp(type, 'p')
            px = dpx;
            py = dpy;
            pz = dpz;
        elseif strcmp (type, 's')
            px = dsx;
            py = dsy;
            pz = dsz;
        else
            error('Wrong type')
        end
        % Calculate transmission coefficients (based on velocity ratios)
        % for P- and S-waves and choose appropriate coefficient for the
        % specified wave type.
        [tpp, tps] = T_fluid_solid(d1, cp1, d2, cp2, cs2, ang1);
        if strcmp(type, 'p')
            T = tpp;
        elseif strcmp(type, 's')
            T = tps;
        end
        % Form up the directivity term
        argx = k1.*uxt.*dx/2;
        argx = argx + eps.*(argx == 0);
        argy = k1.*uyt.*dy/2;
        argy = argy + eps.*(argy == 0);
        dir = (sin(argx)./argx).*(sin(argy)./argy);
        % Form up the denominator term
        D1 = r1 + r2.*(c2/c1).*(cosd(ang1)./cosd(ang2)).^2;
        D2 = r1 + r2.*(c2/c1);
        % Put transmission coefficient, polarization, directivity, phase
        % term and denominator together to calculate velocity components.
        vx = vx + T.*px.*dir.*exp(1i.*k1.*r1 + 1i.*k2.*r2)./sqrt(D1.*D2);
        vy = vy + T.*py.*dir.*exp(1i.*k1.*r1 + 1i.*k2.*r2)./sqrt(D1.*D2);
        vz = vz + T.*pz.*dir.*exp(1i.*k1.*r1 + 1i.*k2.*r2)./sqrt(D1.*D2);
    end
end
% Include external factor for these components
vx = vx.*(-1i*k1*dx*dy)/(2*pi);
vy = vy.*(-1i*k1*dx*dy)/(2*pi);
vz = vz.*(-1i*k1*dx*dy)/(2*pi);

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.