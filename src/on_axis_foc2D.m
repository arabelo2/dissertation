function p = on_axis_foc2D(b, R, f, c, z)
% Ensure no division by zero at z = 0
z = z + eps*(z == 0);

% Define transducer wave number
kb = 2*pi*f*b/c;
% Define u and prevent division by zero
u = (1 - z/R);
u = u + eps*(u == 0);
% Argument of the Fresnel integral and denominator in on-axis pressure
% equation.
x = sqrt((u.*kb.*b)./(pi.*z)).*(z <= R) + sqrt((-u.*kb.*b)./(pi.*z)).*(z > R);
denom = sqrt(u).*(z <= R) + sqrt(-u).*(z > R);
Fr = fresnel_int(x).*(z <= R) + conj(fresnel_int(x)).*(z > R);
% Compute normalized on-axis pressure (p/rho*c*v0) with the propagation
% phase term exp(ikz) removed. Use analytical values near the focus and the
% numerical Fresnel integral values away from the focus.
p = (sqrt(2/1i).*sqrt((b/R).*kb/pi)).*(abs(u) <= .005) + (sqrt(2/1i).*Fr./denom).*(abs(u) > .005);

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.