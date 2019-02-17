function amp = discrete_windows(M, type)
m = 1:M;
switch type
    case 'cos'
        amp = sin(pi*(m - 1)/(M - 1));
    case 'Han'
        amp = (sin(pi*(m - 1)/(M - 1))).^2;
    case 'Ham'
        amp = 0.54 - 0.46*cos(2*pi*(m - 1)/(M -1));
    case 'Blk'
        amp = 0.42 - 0.50*cos(2*pi*(m - 1)/(M - 1)) + 0.08*cos(4*pi*(m - 1)/(M - 1));
    case 'tri'
        amp = 1 - abs(2*(m - 1)/(M - 1) - 1);
    case 'rect'
        amp = ones(1, M);
    otherwise
        disp('Wrong type. Choices are ''cos'', ''Han'', ''Ham'', ''Blk'', ''tri'', ''rect''')
end