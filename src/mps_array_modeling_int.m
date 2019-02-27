% This script solves for the normalized pressure wave field of a 2-D array
% of rectangular elements radiating waves in a fluid using the MATLAB
% function ps_3Dv. Both time delay and apodization laws can specified for
% the array to steer it and focus it.
clear all;
format longG;
% --- Give input parameters ---
tic     % Time the calculation
lx = .15e-3;    % Element length in x-direction
ly = .15e-3;    % Element length in y-direction
gx = .05e-3;    % Gap length in x-direction
gy = .05e-3;    % Gap length in y-direction
f = 5e6;    % Frequency (Hz)
d1 = 1.0;   % Density, medium one (arbitrary units)
cp1 = 1480; % Compressional wave speed, medium one (m/s)
d2 = 7.9;   % Density, medium two (arbitrary units)
cp2 = 5900; % Compressional wave speed, medium two (m/s)
cs2 = 3200; % Shear wave speed, medium two (m/s)
type = 'p'; % Wave type, medium two
mat = {d1, cp1, d2, cp2, cs2, type};    % Form material vector
L1 = 32;    % Number of elements in x-direction
L2 = 1; % Nuber of elements in y-direction
angt = 10.217;  % Angle of the array (deg)
Dt0 = 50.8e-3; % Height of array center from interface (m)
theta2 = 0;  % Steering angle in theta direction (deg)
phi = 0;    % Steering angle in phi direction (deg)
DF = inf;   % Focal distance (m)
lambda = mat{2}/f;
% Weighting choices are 'rect', 'cos', 'Han', 'Ham', 'Blk', 'tri'.
ampx_type = 'rect'; % Weighting coefficients in x-direction
ampy_type = 'rect'; % Weighting coefficients in y-direction

% Field points (x, y, z) to evaluate
N = 10;
xmin = -30e-3;
xmax = 30e-3;
xnpoints = N*ceil(abs(xmax - xmin)/lambda);
xs = linspace(xmin, xmax, xnpoints);

zmin = 1e-3;
zmax = 100e-3;
znpoints = N*ceil(abs(zmax - zmin)/lambda);
zs = linspace(zmin, zmax, znpoints);

y = 0;
[x, z] = meshgrid(xs, zs);
% --- End input parameters ---
c1 = cp1;
if strcmp(type, 'p')
    c2 = cp2;
elseif strcmp(type, 's')
    c2 = cs2;
else
    error ('Type incorrect')
end
% Calculate array pitches
sx = lx + gx;
sy = ly + gy;
% Compute centroid locations for the elements
Nx = 1:L1;
Ny = 1:L2;
ex = (2*Nx - 1 - L1)*sx/2;
ey = (2*Ny - 1 - L2)*sy/2;
% Generate time delays, put in exponential and calculate amplitude weights
td = delay_laws3D_int(L1, L2, sx, sy, angt, phi, theta2, Dt0, DF, c1, c2, 'y');
delay = exp(1i.*2.*pi.*f.*td);
Cx = discrete_windows(L1, ampx_type);
Cy = discrete_windows(L2, ampy_type);
% Calculate normalized velocity
vx = 0;
vy = 0;
vz = 0;
for nn = 1:L1
    for ll = 1:L2
        [vxe, vye, vze] = ps_3Dint(lx, ly, f, mat, ex(nn), ey(ll), angt, Dt0, x, y, z, 1, 1);
        vx = vx + Cx(nn)*Cy(ll)*delay(nn, ll)*vxe;
        vy = vy + Cx(nn)*Cy(ll)*delay(nn, ll)*vye;
        vz = vz + Cx(nn)*Cy(ll)*delay(nn, ll)*vze;
    end
end
%--- Output ---
% Plot results
vmag = sqrt(abs(vx).^2 + abs(vy).^2 + abs(vz).^2);
imagesc(xs, zs, vmag)
shading interp
colormap(jet)
colorbar
axis vis3d
toc % End of time calculations

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.