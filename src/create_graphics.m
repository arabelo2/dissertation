max_Ppp_r = max(max(Ppp_r));
max_Ppp_c = max(max(Ppp_c));
wait = .3;
close all
hold on
for index = (300:5:400)/1000
   zp = length(z(z <= index));
   plot(x, Ppp_r(:, zp)/ max_Ppp_r, '+')
   pause(wait)
   plot(x, Ppp_r(:, zp)/ max_Ppp_r)
   pause(wait)
   plot(x, Ppp_c(:, zp)/ max_Ppp_c, 'o')
   pause(wait)
   plot(x, Ppp_c(:, zp)/ max_Ppp_c)
   pause(wait)
end
hold off

wait = 5;
for index = (0)/1000
   xp = length(x(x <= index));
   %plot(z, Ppp_r(xp, :)/ max_Ppp_r, '+')
   %pause(wait)
   plot(z*1000, Ppp_r(xp, :)/ max_Ppp_r)
   pause(wait)
   %plot(z, Ppp_c(xp, :)/ max_Ppp_c, 'o')
   %pause(wait)
   plot(z*1000, Ppp_c(xp, :)/ max_Ppp_c)
   pause(wait)
end
hold off
grid on
grid minor

a_r=diff(Ppp_r/max_Ppp_r);
a_c=diff(Ppp_c/max_Ppp_c);
dist=norm(a_r - a_c);
tol=0.1;
if dist<tol
   disp('graphs have similar charact.');
else
   disp('graphs have differente charact.');
end