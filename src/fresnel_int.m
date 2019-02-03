function y = fresnel_int(x)
% Separate arguments into positive and negative values, change sign of the
% negative values
xn = -x(x < 0);
xp = x(x >= 0);
% Compute cosine and sine integrals of the negative values, using the
% oddness property of the cosine and sign integrals
[cn, sn] = cs_int(xn);
cn = -cn;
sn = -sn;
% Compute consine and sing integrals of the positive values
[cp, sp] = cs_int(xp);
% Combine cosine and sine integrals for positive and negative values and
% return the complex Fresnel integral
ct = [cn cp];
st = [sn sp];
y = ct + 1i*st;

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.

