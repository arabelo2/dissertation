% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Region: On Edge (Beam: y = R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function omega = OnEdge (c, x, R, t0, t, t1, t2)

% Version 2.0
omega = zeros(1, length(t));
index3 = find(t0 == t & t == t1);
index4 = find(t1 < t & t <= t2);

% if t = t0 = t1
omega(index3) = pi*ones(1, length(index3)); 

% if t1 < t <= t2
omega(index4) = 2*acos(sqrt(c^2*t(index4).^2 - x^2)/(2*R));

% Version 1.0
% n = 1;
% nmax = length(t);
% omega = zeros(1,nmax);
% while ((n <= nmax) && (t(n) < t0))
%     omega(n) = 0;
%     n = n + 1;
% end
% while ((n <= nmax) && (t0 == t(n)) && (t(n) == t1))
%     omega(n) = pi;
%     n = n + 1;
% end
% while ((n <= nmax) && (t1 < t(n)) && (t(n) <= t2))
%     omega(n) = 2*acos(sqrt(c^2*t(n).^2 - x^2)/(2*R));
%     n = n + 1;
% end
% while (n <= nmax)
%     omega(n) = 0;
%     n = n + 1;
% end