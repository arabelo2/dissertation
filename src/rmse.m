max_Ppp_r = max(max(Ppp_r));
max_Ppp_c = max(max(Ppp_c));

for xp = 1:size(Ppp_r, 1)
    rmse_z(xp) = rms(Ppp_c(xp, :)/ max_Ppp_c - Ppp_r(xp, :)/ max_Ppp_r);
end

for zp = 1:size(Ppp_r, 2)
    rmse_x(zp) = rms(Ppp_c(:, zp)/ max_Ppp_c - Ppp_r(:, zp)/ max_Ppp_r);
end

figure()
hold on
xp = length(x(x <= 0/1000))
plot(Ppp_r(xp, :)/ max_Ppp_r)
plot(Ppp_c(xp, :)/ max_Ppp_c)
plot(10*rmse_x./(Ppp_c(xp, :)/max_Ppp_c))
plot(10*rmse_x./(Ppp_r(xp, :)/max_Ppp_r))
hold off

figure()
hold on
zp = length(z(z <= 380/1000))
plot(Ppp_r(:, zp)/ max_Ppp_r)
plot(Ppp_c(:, zp)/ max_Ppp_c)
plot(rmse_z./(Ppp_r(:, zp)'/max_Ppp_r)); hold on
plot(rmse_z./(Ppp_c(:, zp)'/max_Ppp_c))
hold off

% http://kawahara.ca/root-mean-square-error-tutorial-matlab/