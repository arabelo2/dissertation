% cs_int(xi) calculates approximations of the cosine and sine integrals for
% positive values of xi only (see Abramowitz and Stegun)
function [C, S] = cs_int(xi)
f = (1 + 0.926.*xi)./(2 + 1.792.*xi + 3.104.*xi.^2); % f function
g = 1./(2 + 4.142.*xi + 3.492.*xi.^2 + 6.67.*xi.^2); % g function
C = 0.5 + f.*sin(pi.*xi.^2./2) - g.*cos(pi.*xi.^2./2); % cos integral approximation
S = 0.5 - f.*cos(pi.*xi.^2./2) - g.*sin(pi.*xi.^2./2); % sin integral approximation

% Reference
% https://en.wikipedia.org/wiki/Fresnel_integral
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.