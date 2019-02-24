function xi = pts_3Dint(ex, ey, xn, yn, angt, Dt0, c1, c2, x, y, z)
% Calculate wave speed ratio
cr = c1/c2;
% Determine size of array needed for xi calculations based on the sizes of
% the (x, y, z) variables and also determine those sizes.
[xi, P, Q] = init_xi3D(x, y, z);
[nrx, ncx] = size(x);
[nry, ncy] = size(y);
[nrz, ncz] = size(z);
% Call ferrari2 function to compute xi with the arguments of that function
% determined by the sizes of the (x, y, z) variables.
De = Dt0 + (ex + xn)*sind(angt);
for pp = 1:P
    for qq = 1:Q
        % x and y are points, z is a row or column vector
        if nrx == 1 && ncx == 1 && nry == 1 && ncr == 1
            Db = sqrt((x - (ex + xn)*cosd(angt)).^2 + (y - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z(pp, qq), De, Db);
        % y and z are points, x is a row or column vector    
        elseif nry == 1 && ncy == 1 && nrz == 1 && ncz == 1
            Db = sqrt((x(pp, qq) - (ex + xn)*cosd(angt)).^2 + (y  - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z, De, Db);
        % x and z are points, y is a row or column vector
        elseif nrx == 1 && ncx == 1 && nrz == 1 && ncz == 1
            Db = sqrt((x - (ex + xn)*cosd(angt)).^2 + (y(pp, qq) - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z, De, Db);
        % y is a point, x and z are equal size PxQ matrices
        elseif nry == 1 && ncy == 1 && nrx == nrz && ncx == ncz
            Db = sqrt((x(pp, qq) - (ex + xn)*cosd(angt)).^2 + (y - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z(pp, qq), De, Db);
        % z is a point, x and y are equal size PxQ matrices
        elseif nrz == 1 && ncz == 1 && nrx == nry && ncx == ncy
            Db = sqrt((x(pp, qq) - (ex + xn)*cosd(angt)).^2 + (y(pp, qq) - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z, De, Db);
        % x is a point, y and z are equal size PxQ matrices
        elseif nrx == 1 && ncx == 1 && nry == nrz && ncy == ncz
            Db = sqrt((x - (ex + xn)*cosd(angt)).^2 + (y(pp, qq) - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z(pp, qq), De, Db);
        % x, y and z are all equal size row and column vectors
        else
            Db = sqrt((x(pp, qq) - (ex + xn)*cosd(angt)).^2 + (y(pp, qq) - (ey + yn)).^2);
            xi(pp, qq) = ferrari2(cr, z(pp, qq), De, Db); 
        end
    end    
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.