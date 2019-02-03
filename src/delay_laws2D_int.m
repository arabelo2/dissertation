function td = delay_laws2D_int(M, s, angt, ang20, DT0, DF, c1, c2, plt)
cr = c1/c2; % Wave speed ratio
Mb = (M -1)/2;
% Compute locaton of element centroids, e;
m = 1:M;
e = (m - 1 - Mb)*s;
% Computed parameters:
% ang10 - incident angle of central ray, deg;
% DXO - distance along interface from center of array to focal point;
% DT - heights of elements above interface;
% DX - distances along interface from elements to focal point;
ang10 = asind(cr.*sind(ang20));
DX0 = DF.*tand(ang20) + DT0.*tand(ang10);
DT = DT0 + e.*sind(angt);
DX = DX0 - e.*cosd(angt);
switch(DF)
    % Steering only case, use linear law.
    case inf
        if (ang10 - angt) > 0
            td = (m - 1)*s*sind(ang10 - angt)/c1;
        else
            td = (M - m)*s*abs(sind(ang10 - angt))/c1;
        end
        % Plotting rays option
        if strcmp(plt, 'y')
            for nn = 1:M
                xp2(1, nn) = e(nn)*cosd(angt);
                yp2(1, nn) = DT(nn);
                xp2(2, nn) = e(nn)*cosd(ang10 - angt)/cosd(ang10) + DT0*tand(ang10);
                dm = e(nn)*cosd(ang10 - angt)/cosd(ang10);
                if ang20 > 0
                    dM = e(M)*cosd(ang10 - angt)/cosd(ang10);
                else
                    dM = e(1)*cosd(ang10 - angt)/cosd(ang10);
                end
                yp2(2, nn) = 0;
                xp2(3, nn) = xp2(2, nn) + (dM - dm)*sind(ang20)*sind(ang20);
                yp2(3, nn) = -(dM -dm)*sind(ang20)*cosd(ang20);
            end
            plot(xp2, yp2, 'b')
        end
        % End plotting rays option
        % Steering and focusing case
    otherwise
        % Solve for ray intersection locations on interface, xi, and path
        % lengths in medium 1 and medium 2, r1 and r2.
        xi = zeros(1, M);
        r1 = zeros(1, M);
        r2 = zeros(1, M);
        for mm = 1:M
            xi(mm) = ferrari2(cr, DF, DT(mm), DX(mm));
            r1(mm) = sqrt(xi(mm)^2 + (DT0 + e(mm)*sind(angt))^2);
            r2(mm) = sqrt((xi(mm) + e(mm)*cosd(angt) - DX0)^2 + DF^2);
        end
        % Solve for time advances, turn into delays, and make the delays,
        % td, positive.
        t = r1/c1 + r2/c2;
        td = max(t) - t
        % Plotting rays option
        if strcmp(plt, 'y')
            for nn = 1:M
                xp(1, nn) = e(nn)*cos(angt*pi/180);
                yp(1, nn) = DT(nn);
                xp(2, nn) = e(nn)*cos(angt*pi/180) + xi(nn);
                yp(2, nn) = 0;
                xp(3, nn) = DX0;
                yp(3, nn) = -DF;                
            end
            plot(xp, yp, 'b')
        end
        % End plotting rays option
end
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.