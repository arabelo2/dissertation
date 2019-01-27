f0 = 5e6
c = 1500
b = 2e-3/2
e = 4e-4
x = (-40:0.1:40)*1e-3;
z = (40)*1e-3

lambda = c/f0

if 2*b > lambda/10
    Nopt = ceil(20*f0*b/c);
else
    Nopt = 1;
end

Nopt

% p = rs_2Dv(b, f, c, e, x, z, Nopt);
% plot(z, abs(p))
% hold on
% plot(z, abs(p), 'ko')
% hold off


p = rs_2Dv(b, f0, c, e, x, z, Nopt);
plot(rad2deg(asin(x/z)), abs(p))
hold on
plot(rad2deg(asin(x/z)), abs(p), 'ko')
hold off

% p = cell(size(x));
% for xx = 1:length(x)
%     p{xx} = rs_2Dv(b, f, c, e, x(xx), z, Nopt);
% end
% 
% for xx = 1:length(x)
%     plot(z, abs(p{xx}))
%     pause(.5)
% end