function td = delay_laws3D(M, N, sx, sy, theta, phi, F, c)
m = 1:M;
n = 1:N;
Mb = (M - 1)/2;
Nb = (N - 1)/2;
exm = (m - 1 - Mb)*sx;
eyn = (n - 1 - Nb)*sy;

% Calculate delays
switch(F)
    % If steering only specified, use explicit steering law
    case(inf)
        for mm = 1:M
            for nn = 1:N
                dt(mm, nn) = (exm(mm)*sind(theta)*cosd(phi) + eyn(nn)*sind(theta)*sind(phi))/c;
            end
        end
        % Make delays all positive
        td = abs(min(min(dt))) + dt;
    % Otherwise, if steering and focusing specified, use time delays to the specified point
    otherwise
        for mm = 1:M
            for nn = 1:N
                r(mm, nn) = sqrt((F*sind(theta)*cosd(phi) - exm(mm))^2 + (F*sind(theta)*sind(phi) - eyn(nn))^2 + F^2*(cosd(theta))^2);
            end
        end
        td = max(max(r/c)) - r/c;
end
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.