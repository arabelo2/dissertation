function p = fresnel_2D (b, f, c, x, z)
% Calculate wave number
kb = 2*pi*f*b/c;
% Put (x, z) coordinates in normalized form
xb = x/b;
zb = z/b;
% Calculate term in Fresnel integral argument
arg = sqrt(kb./(pi*zb));
% Calculate normalized pressure
p = sqrt(1/(2*1i)).*exp(1i*kb*zb).*(fresnel_int(arg.*(xb + 1)) - fresnel_int(arg.*(xb -1)));