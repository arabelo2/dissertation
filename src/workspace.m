f = 5e6
c = 1500
b = 6.35/2
e = b/100
x = 0
z = 0:100

lambda = c/f;

if 2*b*1e-3 > lambda/10
    Nopt = ceil(20e-3*f*b/c);
else
    Nopt = 1
end

Nopt
2e-3*b/lambda

p = rs_2Dv(b, f, c, e, x, z, Nopt);
plot(z, abs(p))
hold on
plot(z, abs(p), 'ro')
