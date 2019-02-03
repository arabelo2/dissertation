function p = Gauss_2D(b, f, c, x, z)
% Retrieve Wen and Breazeale coefficients
[A ,B] = gauss_c15;
% Calculate the wave number
kb = 2*pi*f*b/c;
% Normalize the (x, z) coordinates
xb = x/b;
zb = z/b;
% Initialize the pressure to zero and then superimpose 15 Gaussian beams to
% calculate the pressure wave field.
p = 0;
for nn = 1:15
    qb = zb - 1i*pi*f*b./(B(nn)*c);
    qb0 = -1i*pi*f*b./(B(nn)*c);
    p = p + sqrt(qb0./qb).*A(nn).*exp(1i*kb*xb.^2./(2*qb));
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.