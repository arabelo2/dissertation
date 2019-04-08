function [A, d, g, xc] = elements(f, c, dl, gd, N)
% dl is the element diameter, d, divided by the wavelength, l, i.e. dl = d/l.
d = dl.*c./f;
% gd is the gap size, g, between elements as a fraction of the element
% size, i.e. gd = g/d.
g = gd.*d;
% A is the total aperture size of the array.
A = N*d + (N - 1)*g;
% x = xc is the location of the centroid of each element where x = 0 is at
% the center of the array.
for nn = 1:N
    xc(nn) = (g + d)*((2*nn - 1)/2 - N/2);
end