function omega = OnEdge (c, x, R, t0, t, t1, t2)
n = 1;
nmax = length(t);
omega = zeros(1,nmax);
while ((n <= nmax) && (t(n) < t0))
    omega(n) = 0;
    n = n + 1;
end
while ((n <= nmax) && (t0 == t(n)) && (t(n) == t1))
    omega(n) = pi;
    n = n + 1;
end
while ((n <= nmax) && (t1 < t(n)) && (t(n) <= t2))
    omega(n) = 2*acos(sqrt(c^2*t(n).^2 - x^2)/(2*R));
    n = n + 1;
end
while (n <= nmax)
    omega(n) = 0;
    n = n + 1;
end