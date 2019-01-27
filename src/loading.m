tic()
f0 = 5e6;
c = 1480;
width = 6e-3/2;
elements_x = 1;
height = 2*width;
elements_y = 1;
kerf = 4e-4;
lambda = c/f0;
xmin = -15e-3; %-(3*width + kerf) * (elements_x/2+1);
xmax = +15e-3; %(3*width + kerf) * (elements_x/2+1);
ymin = 0;
ymax = 0;
zmin = 0;
zmax = +80e-3; % 400 * lambda;
focus_x = 0;
focus_y = 0;
focus_z = 25 * lambda;
xpoints = 512;
ypoints = 1;
zpoints = 512;
dx = (xmax-xmin)/xpoints;
dy = (ymax-ymin)/ypoints;
dz = (zmax-zmin)/zpoints;
x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;

% x = 0;
% z  = 60e-3;

if 2*width > lambda/10
    Nopt = ceil(20*f0*width/c);
else
    Nopt = 1;
end

P1 = cell(length(z), length(x));
P2 = cell(length(z), length(x));

% % for xx = 1:length(x)
% %     for zz = 1:length(z)
% %         P1{zz, xx} = rs_2Dv(width, f0, c, kerf, x(xx), z(zz), Nopt);
% %     end
% % end

% p = rs_2Dv(width, f0, c, kerf, x, z, Nopt);
% plot(rad2deg(asin(x./z)), abs(p))
%%%%% plot(z, abs(p))

% % figure(1)
% % pcolor(x, z, cellfun(@abs, P))
% % shading interp
% % colormap(jet)
% % axis vis3d
% % grid on
% % grid minor
% % set(gca,'FontSize',20);

% % figure(2)
% % mesh(x, z, cellfun(@abs, P))
% % shading interp
% % colormap(jet)
% % axis normal
% % grid on
% % grid minor
% % set(gca,'FontSize',20);


mat = [1000, 1500, 1000, 1500];
angt = 0;
Dt0 = 1e-3;
p = ls_2Dint(width, f0, mat, kerf, angt, Dt0, x, z, Nopt);

for xx = 1:length(x)
    for zz = 1:length(z)
        P1{zz, xx} = rs_2Dv(width, f0, c, kerf, x(xx), z(zz), Nopt);
        P2{zz, xx} = ls_2Dint(width, f0, mat, kerf, angt, Dt0, x, z, Nopt);
    end
end
