% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Region: Inside geometrical (Beam: y < R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function omega = InsideGeometrical (c, x, y, R, t0, t, t1, t2)

% Version 2.0
omega = zeros(1, length(t));
index1 = find(t0 <= t &  t <= t1);
index2 = find(t1 < t  &  t <= t2);

% if t0 <= t <= t1
omega(index1) = 2*pi*ones(1, length(index1)); 

% if t1 < t <= t2
omega(index2) = 2*acos((c^2*t(index2).^2 - x^2 + y^2 - R^2)./(2*y*sqrt(c^2*t(index2).^2 - x^2))); 

% Version 1.0
% n = 1;
% nmax = length(t);
% omega = zeros(1, nmax);
% while ((n <= nmax) && (t(n) < t0))
%     omega(n) = 0;
%     n = n + 1;
% end
% while ((n <= nmax) && (t0 <= t(n)) && (t(n) <= t1))
%     omega(n) = 2*pi;
%     n = n + 1;
% end
% while ((n <= nmax) && (t1 < t(n)) && (t(n) <= t2))
%     omega(n) = 2*acos((c^2*t(n).^2 - x^2 + y^2 - R^2)./(2*y*sqrt(c^2*t(n).^2 - x^2)));
%     n = n + 1;
% end
% while (n <= nmax)
%     omega(n) = 0;
%     n = n + 1;
% end