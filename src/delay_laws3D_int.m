function td = delay_laws3D_int(Mx, My, sx, sy, thetat, phi, theta2, DT0, DF, c1, c2, plt)
% Compute wave speed ratio
cr = c1/c2;
% Compute element centroid locations
mx = 1:Mx;
my = 1:My;
Mbx = (Mx - 1)/2;
Mby = (My - 1)/2;
ex = (mx - 1 - Mbx)*sx;
ey = (my - 1 - Mby)*sy;
% Initialize variables to be used
t = zeros(Mx, My);
Db = zeros(Mx, My);
De = zeros(1, Mx);
xi = zeros(Mx, My);
% Angle in first medium (Deg)
ang1 = asind(c1*sind(theta2)/c2);
% Calculate delays
switch(DF)
    % If steering only specified, use explicit steering law
    case(inf)
        ux = sind(ang1)*cosd(phi)*cosd(thetat) - cosd(ang1)*sind(thetat);
        uy = sind(ang1)*sind(phi);
        for m = 1:Mx
            for n = 1:My
                t(m, n) = (ux*ex(m) + uy*ey(n))/c1;
            end
        end
        % Make delays all positive
        td = abs(min(min(t))) + t;
    % Otherwise, if steering and focusing specified, use time delays to the specified point
    otherwise
        % Determine distances De and Db needed in arguments of ferrari2
        % function.
        DQ = DT0*tand(ang1) + DF*tand(theta2);
        x = DQ*cosd(phi);
        y = DQ*sind(phi);
        for m = 1:Mx
            for n = 1:My
                Db(m, n) = sqrt((x - ex(m)*cosd(thetat)).^2 + (y -ey(n)).^2);
            end
        end
        De = DT0 + ex*sind(thetat);
        % Use ferrari2 method to determine distance, xi, where a ray from
        % an element to the point (x, y, DF) intersects the interface in
        % the plane of incidence.
        for m = 1:Mx
            for n = 1:My
                xi(m, n) = ferrari2(cr, DF, De(m), Db(m, n));
            end
        end
        % Use ray distances to calculate time advances
        for m = 1:Mx
            for n = 1:My
                t(m, n) = sqrt(xi(m, n).^2 + De(m).^2)/c1 + sqrt(DF^2 + (Db(m, n) - xi(m, n)).^2)/c2;
            end
        end        
        % Turn time advances into delays and make all delays positive
        td = max(max(t)) - t;
        % Plotting rays option
        if strcmp(plt, 'y')
            for m = 1:Mx
                for n = 1:My
                xp(1, 1) = ex(m)*cosd(thetat);
                zp(1, 1) = DT0 + ex(m)*sind(thetat);
                yp(1, 1) = ey(n);
                xp(2, 1) = ex(m)*cosd(thetat) + xi(m, n)*(x - ex(m)*cosd(thetat))/Db(m, n);
                yp(2, 1) = ey(n) + xi(m, n)*(y - ey(n))/Db(m, n);
                zp(2, 1) = 0;
                xp(3, 1) = x;
                yp(3, 1) = y;
                zp(3, 1) = -DF;
                plot3(xp, yp, zp)
                [xp yp zp]
                pause(5)
                hold on
                end
            end
            hold off
        end % End plotting rays option
end
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.