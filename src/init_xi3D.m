function [xi, P, Q] = init_xi3D(x, y, z)
% Get sizes of (x, y, z)
[nrx, ncx] = size(x);
[nry, ncy] = size(y);
[nrz, ncz] = size(z);
% If x and z are equal sized matrices [nrx, ncx] and y is single value,
% make xi a [nrx, ncx] matrix.
if nrx == nrz && ncx == ncz && nry == 1 && ncy == 1
    xi = zeros (nrx, ncx);
    P = nrx;
    Q = ncx;    
% If x and y are equal sized matrices [nrx, ncx] and z is single value,
% make xi a [nrx, ncx] matrix.
elseif nrx == nry && ncx == ncy && nrz == 1 && ncz == 1
    xi = zeros (nrx, ncx);
    P = nrx;
    Q = ncx;
% If y and z are equal sized matrices [nry, ncy] and x is single value,
% make xi a [nry, ncy] matrix.
elseif nry == nrz && ncy == ncz && nrx == 1 && ncx == 1
    xi = zeros (nry, ncy);
    P = nry;
    Q = ncy;
% If z is a [1, ncz] vector and x and y are single value, make xi a [1, ncz] vector
elseif nrz == 1 && ncz > 1 && nrx == 1 && ncx == 1 && nry == 1 && ncy == 1
    xi = zeros(1, ncz);
    P = 1;
    Q = ncz;
% If z is a [nrz, 1] vector and x and y are single values, make xi a [nrz, 1] vector
elseif ncz == 1 && nrz > 1 && nrx == 1 && ncx == 1 && nry == 1 && ncy == 1
    xi = zeros(nrz, 1);
    P = nrz;
    Q = 1;
% If x is a [1, ncx] vector and y and z are single values, make xi a [1, ncx] vector
elseif nrx == 1 && ncx > 1 && nry == 1 && ncy == 1 && nrz == 1 && ncz == 1
    xi = zeros(1, ncx);
    P = 1;
    Q = ncx;
% If x is a [nrx, 1] vector and y and z are single values, make xi a [nrx, 1] vector
elseif ncx == 1 && nrx > 1 && nry == 1 && ncy == 1 && nrz == 1 && ncz == 1
    xi = zeros(nrx, 1);
    P = nrx;
    Q = 1;
% If y is a [1, ncy] vector and x and z are single values, make xi a [1, ncy] vector
elseif nry == 1 && ncy > 1 && nrx == 1 && ncx == 1 && nrz == 1 && ncz == 1
    xi = zeros(1, ncy);
    P = 1;
    Q = ncy;
% If y is a [nry, 1] vector and x and z are single values, make xi a [nry, 1] vector
elseif ncy == 1 && nry > 1 && nrx == 1 && ncx == 1 && nrz == 1 && ncz == 1
    xi = zeros(nry, 1);
    P = nry;
    Q = 1;
% If x, y and z are equal size [1, ncx] vectors, make xi a [1, ncx] vector
elseif nrx == nry && ncx == ncy && nrz == nrx && ncz == ncx && nrx == 1
    xi = zeros(1, ncx);
    P = 1;
    Q = ncx;
% If x, y and z are equal size [nrx, 1] vectors, make xi a [nrx, 1] vector
elseif nrx == nry && ncx == ncy && nrz == nrx && ncz == ncx && ncx == 1
    xi = zeros(nrx, 1);
    P = nrx;
    Q = 1;
% Other combinations are not supported
else
    error('(x, y, z) combination given is not supported.')
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.