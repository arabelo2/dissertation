function [tpp, tps] = T_fluid_solid(d1, cp1, d2, cp2, cs2, theta1)
% Put incident angle in radians
iang = (theta1.*pi)./180;
% Calculate sin(theta) for refracted p- and s-waves
sinp = cp2/cp1*sin(iang);
sins = cs2/cp1*sin(iang);
% Calculate cos(theta) for refracted p- and s-waves for angles beyond
% critical, the value of the cosine is computed for positive frequencies
% only.
cosp = 1i*sqrt(sinp.^2 - 1).*(sinp >= 1) + sqrt(1 - sinp.^2).*(sinp < 1);
coss = 1i*sqrt(sins.^2 - 1).*(sins >= 1) + sqrt(1 - sins.^2).*(sins < 1);
% Calculate transmission coefficients
denom = cosp + d2/d1*cp2/cp1*sqrt(1 - sin(iang).^2).*(4.*((cs2/cp2)^2).*(sins.*coss.*sinp.*cosp) + 1 -4*(sins.^2).*(coss.^2));
tpp = (2*sqrt(1 - sin(iang).^2).*(1 - 2*(sins.^2)))./denom;
tps = -(4*cosp.*sins.*sqrt(1 - sin(iang).^2))./denom;

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.