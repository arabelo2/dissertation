% Transforms corresponding elements of the focal point spherical coordinate arrays azimuth (theta), elevation (psi), and F to Cartesian.

% \PHII     --> Rotate around the z-axis (Roll)
% \THETA    --> Rotate around the y-axis (Pitch) / Azimuth is the counterclockwise angle in the z-x plane measured in radians from the positive z-axis.
% \PSI      --> Rotate around the x-axis (Yaw) / Elevation is the elevation angle in radians from the z-x plane.
% F - Focal distance [m]

function [zf, xf, yf] = focalPoint (theta, psi, F)
    [zf, xf, yf] = sph2cart(deg2rad(theta), deg2rad(psi), F);
end