function [xi, P, Q] = init_xi (x, z)
% Get sizes of x and z variables
[nrx, ncx] = size(x);
[nrz, ncz] = size(z);
% If x and z are equal sized matrices, vectors, or scalars, xi is of the
% same size
if nrx == nrz && ncx == ncz
    xi = zeros (nrx, ncx);
    P = nrx;
    Q = ncx;
% If x is a column vector and z a scalar, xi is the same size column vector
elseif nrx > 1 && ncx == 1 && nrz == 1 && ncz == 1
    xi = zeros(nrx, 1);
    P = nrx;
    Q = 1;
% If z is a column vector and x a scalar, xi is the same size column vector
elseif nrz > 1 && ncz == 1 && nrx == 1 && ncx == 1
    xi = zeros(nrz, 1);
    P = nrz;
    Q = 1;
% If x is a row vector and z a scalar, xi is the same size row vector
elseif ncx > 1 && nrx == 1 && nrz == 1 && ncz == 1
    xi = zeros(1, ncx);
    P = 1;
    Q = ncx;
% If z is a row vector and x a scalar, xi is the same size row vector
elseif ncz > 1 && nrz == 1 && nrx == 1 && ncx == 1
    xi = zeros(1, ncz);
    P = 1;
    Q = ncz;
% Other combinations are not supported
else
    error('(x, z) must be (vector, scalar) pairs or equal matrices')
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.