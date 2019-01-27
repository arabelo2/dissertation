function xi = pts_2Dintf(e, xn, angt, Dt0, c1, c2, x, z)
% Calculate wave speed ratio
cr = c1/c2;
% Based on sizes of (x, z), determine corresponding number of rows and
% columns (P, Q) needed for xi calculations and initialize xi as zeros.
[xi, P, Q] = init_xi(x, z);
% Obtain sizes of (x, z) so appropriate arguments can be found in the calls
% to the function ferrari2 when making the xi calculations
[nrx, ncx] = size (x);
[nrz, ncz] = size (z);
% Calculate xi locations using Ferrari's method
for pp = 1:P
    for qq = 1:Q
        Dtn = Dt0 + (e + xn)*sin(angt*pi/180);
        % If x is a point, and z is a rwo or column vector
        if nrx == 1 && ncx == 1
            Dxn = x - (e + xn)*cos(angt*pi/180);
            xi(pp, qq) = ferrari2(cr, z(pp, qq), Dtn, Dxn);
        % If z is a point, and x is a row or column vector
        elseif nrz == 1 && ncz == 1
            Dxn = x(pp, qq) - (e + xn)*cos(angt*pi/180);
            xi(pp, qq) = ferrari2(cr, z, Dtn, Dxn);
        % If x and z are equal size PxQ matrices
        else
            Dxn = x(pp, qq) - (e + xn)*cos(angt*pi/180);
            xi(pp, qq) = ferrari2(cr, z(pp, qq), Dtn, Dxn);
        end
    end
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.