max_Ppp = max(max(Ppp'));
max_picoapico = max(max(picoapico));

for xp = 1:size(Ppp, 1)
    rmse_z(xp) = rms(picoapico(xp, :)/ max_picoapico - Ppp(xp, :)'/ max_Ppp);
end

for zp = 1:size(Ppp, 2)
    rmse_x(zp) = rms(picoapico(:, zp)/ max_picoapico - Ppp(:, zp)'/ max_Ppp);
end

figure()
hold on
xp = length(x(x <= 0/1000))
plot(Ppp(xp, :)'/ max_Ppp)
plot(picoapico(xp, :)/ max_picoapico)
plot(10*rmse_x./(picoapico(xp, :)/max_picoapico))
plot(10*rmse_x./(Ppp(xp, :)'/max_Ppp))
hold off

figure()
hold on
zp = length(z(z <= 380/1000))
plot(Ppp(:, zp)'/ max_Ppp)
plot(picoapico(:, zp)/ max_picoapico)
plot(rmse_z./(Ppp(:, zp)'/max_Ppp)); hold on
plot(rmse_z./(picoapico(:, zp)/max_picoapico))
hold off

% http://kawahara.ca/root-mean-square-error-tutorial-matlab/