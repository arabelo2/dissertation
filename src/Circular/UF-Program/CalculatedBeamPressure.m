function [E, P,  Pi, phii, omega, vn, ka, S, t, t0, t1, t2, t_vn, t_pconv, t_econv] = CalculatedBeamPressure(x, c, R, y , delta_t, Uo, f0, rho, T, lambda, ncycle, K);
y = abs(y); % Absolute value
t0 = x/c;
t1 = (1/c)*((R - y)^2 + x^2)^(1/2);
t2 = (1/c)*((R + y)^2 + x^2)^(1/2);
t_vn = 0:delta_t:ncycle*T;

if (t0 - delta_t < 0)
    t = t0:delta_t:t2 + delta_t;
else
    t = t0 - delta_t:delta_t:t2 + delta_t;
end

% The angle of equidistant arc included on the surface.
if (y < R)
    omega = InsideGeometrical (c, x, y, R, t0, t, t1, t2);
elseif (y == R)
    omega = OnEdge (c, x, R, t0, t, t1, t2);
else
    omega = OutsideGeometrical (c, x, y, R, t, t1, t2);
end

% The velocity potential impulse response.
phii = c.*omega/(2*pi);

% The pressure impulse response.
Pi = rho*diff(phii)/delta_t;

% The input velocity of disc surface.
vn = InputVelocityOfDiscSurface (Uo, f0, t_vn);

% The pressure.
P=conv(vn, Pi);

% The outage voltage (The echo response).
E = (-1)*(K*rho/(2*c))*conv(vn, conv(diff(phii)/delta_t, diff(phii)/delta_t));

%Número de Seki
S = x*lambda/R^2;

ka = 2*pi*R/lambda; % Constant.

omega(end:end+numel(t)-numel(omega))=nan;
phii(end:end+numel(t)-numel(phii))=nan;
Pi(end:end+numel(t)-numel(Pi))=nan;

% Time
t_pconv=t0 + delta_t*(1:length(P));
t_econv=2*t0 + delta_t*(1:length(E));