f0 = 5e6;
c1 = 1500;
b = 12e-3/2;
elements_x = 32;
height = 2*b;
elements_y = 1;
e = 4e-4;
lambda = c1/f0;
xmin = -15e-3; %-(3*width + kerf) * (elements_x/2+1);
xmax = +15e-3; %(3*width + kerf) * (elements_x/2+1);
ymin = 0;
ymax = 0;
zmin = 0;
zmax = 100e-3; % 400 * lambda;
focus_x = 0;
focus_y = 0;
focus_z = 100e-3;
xpoints = 1024;
ypoints = 1;
zpoints = 1024;
dx = (xmax-xmin)/xpoints;
dy = (ymax-ymin)/ypoints;
dz = (zmax-zmin)/zpoints;
x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;

x = 0;
% z  = 90e-3;

M = elements_x;
s = 0.5e-3;
angt = 0;
ang20 = 75;
DT0 = 10e-3;
DF = 10e-3;
c2 = 5900;
plt = 'n';
td = delay_laws2D_int(elements_x, s, angt, ang20, DT0, DF, c1, c2, plt)

p = fresnel_2D (b, f0, c1, x, z)
plot(z, abs(p), 'b')
hold on

p = Gauss_2D(b, f0, c1, x, z)
plot(z, abs(p), 'r')
hold on

p = on_axis_foc2D(b, focus_z, f0, c1, z);
plot(z, abs(p))
hold on

% % ylim([0 1.5])

if 2*b > lambda/10
    Nopt = ceil(20*f0*b/c1);
else
    Nopt = 1;
end

P1 = cell(length(z), length(x));
P2 = cell(length(z), length(x));
P3 = cell(length(z), length(x));

% % for xx = 1:length(x)
% %     for zz = 1:length(z)
% %         P1{zz, xx} = rs_2Dv(b, f0, c, e, x(xx), z(zz), Nopt);
% %     end
% % end

% p = rs_2Dv(b, f0, c, e, x, z, Nopt);
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


% % mat = [997, 1480, 997, 1480];
% % angt = 0;
% % Dt0 = 60e-3;
% % p = ls_2Dint(b, f0, mat, e, angt, Dt0, x, z, Nopt);

% for xx = 1:length(x)
%     for zz = 1:length(z)
%         % P1{xx, zz} = rs_2Dv(b, f0, c1, e, x(xx), z(zz), Nopt);
%         % P2{xx, zz} = ls_2Dint(b, f0, mat, e, angt, Dt0, x(xx), z(zz), Nopt);
%         P2{xx, zz} = fresnel_2D (b, f0, c1, x(xx), z(zz));
%     end
% end

% 2D
% figure(3)
% pcolor(z, x, cellfun(@abs, P2))
% shading interp
% colormap(jet)
% axis vis3d
% grid on
% grid minor
% set(gca,'FontSize',20);
% camroll(-90)

% 3D
% figure(4)
% mesh(z, x, cellfun(@abs, P2))
% shading interp
% colormap(jet)
% axis normal
% grid on
% grid minor
% set(gca,'FontSize',20);

