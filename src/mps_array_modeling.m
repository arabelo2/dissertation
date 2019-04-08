% This script solves for the normalized pressure wave field of a 2-D array
% of rectangular elements radiating waves in a fluid using the MATLAB
% function ps_3Dv. Both time delay and apodization laws can specified for
% the array to steer it and focus it.
clear all;
format longG;
% --- Give input parameters ---
lx = .15e-3; % Element length in x-direction
ly = .15e-3; % Element length in y-direction
gx = .05e-3; % Gap length in x-direction
gy = .05e-3; % Gap length in y-direction
f = 5e6; % Frequency (Hz)
c = 1500; % Wave speed (m/s)
L1 = 11; % Number of elements in x-direction
L2 = 11; % Nuber of elements in y-direction
theta = 20; % Steering angle in theta direction (deg)
phi = 0; % Steering angle in phi direction (deg)
F1 = 50.4e-3; % Focal distance (m)
lambda = c/f;
% Weighting choices are 'rect', 'cos', 'Han', 'Ham', 'Blk', 'tri'.
ampx_type = 'rect'; % Weighting coefficients in x-direction
ampy_type = 'rect'; % Weighting coefficients in y-direction

% Field points (x, y, z) to evaluate
N = 10;
xmin = -10e-3;
xmax = +10e-3;
xnpoints = 1000; % xnpoints = N*ceil(abs(xmax - xmin)/lambda);
xs = linspace(xmin, xmax, xnpoints);

% zmin = 0;
% zmax = +120e-3;
% znpoints = 300; % znpoints = N*ceil(abs(zmax - zmin)/lambda);
% zs = linspace(zmin, zmax, znpoints);

ymin = -10e-3;
ymax = +10e-3;
ynpoints = 1000; % ynpoints = N*ceil(abs(ymax - ymin)/lambda);
ys = linspace(ymin, ymax, ynpoints);

z = 4e-3;
% [x, z] = meshgrid(xs, zs);
[x, y] = meshgrid(xs, ys);
% --- End input parameters ---
% Calculate array pitches
sx = lx + gx;
sy = ly + gy;
% Compute centroid locations for the elements
Nx = 1:L1;
Ny = 1:L2;
ex = (2*Nx - 1 - L1)*sx/2;
ey = (2*Ny - 1 - L2)*sy/2;
% Generate time delays, put in exponential and calculate amplitude weights
td = delay_laws3D(L1, L2, sx, sy, theta, phi, F1, c);
delay = exp(1i.*2.*pi.*f.*td);
Cx = 1; % Cx = discrete_windows(L1, ampx_type);
Cy = 1; % Cy = discrete_windows(L2, ampy_type);
% Calculate normalized pressure
p = 0;
for nn = 1:L1
    for mm = 1:L2
        % p = p + Cx(nn)*Cy(mm)*delay(nn, mm)*ps_3Dv(lx, ly, f, c, ex(nn), ey(mm), x, y, z); 
        p = p + Cx*Cy*delay(nn, mm)*ps_3Dv(lx, ly, f, c, ex(nn), ey(mm), x, y, z); 
        imagesc(xs, ys, abs(p))
        shading interp
        colormap(gray) % colormap(jet)
        colorbar
        axis vis3d
        pause(.05)
    end
end
%--- Output ---
% Plot results based on specifications of (x, y, z) points
% imagesc(xs, zs, abs(p))
% shading interp
% colormap(jet)
% colorbar
% axis vis3d

% Plot results based on specifications of (x, y, z) points
imagesc(xs, ys, abs(p))
shading interp
colormap(gray) % colormap(gray)
colorbar
axis vis3d

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.