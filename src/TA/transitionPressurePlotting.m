hold off
figure(50)
Pmmaxt = max(max(real(cellfun(@max, P))));
% figure(10)
for zz = 1:length(z)
% for zz = 40
% for zz = (length(find (z < F*cos(20*pi/180))) - 5):(length(find (z < F*cos(20*pi/180))) + 0)
% for zz = (length(find (z < F)) - 5):(length(find (z < F)) + 5)
    for yy = 1:length(y(1, :))
        for xx = round(length(x(1, :))/ 2)
        % for xx = 1:length(x)
        % for xx = (round(length(x(1, :))/ 2) - 10):(round(length(x(1, :))/ 2) + 10)
            % Ploting individual transient pressure
            plot(P{yy, xx, zz}/ Pmmaxt)
            % plot(P{yy, xx, zz}(find(P{yy, xx, zz}, 1):find(P{yy, xx, zz}, 1, 'last'))/Pmmax, 'ro')
            % plot(P{yy, xx, zz}(find(P{yy, xx, zz} == max(P{yy, xx, zz}))-170: find(P{yy, xx, zz} == max(P{yy, xx, zz})) + 900)/Pmmaxt, 'r')
            % plot(P{yy, xx, zz}(find(P{yy, xx, zz} == max(P{yy, xx, zz})): find(P{yy, xx, zz} == max(P{yy, xx, zz})) + 1024)/Pmmaxt)
            % plot(t_conv{yy, xx, zz}(find(P{yy, xx, zz}, 1):find(P{yy, xx, zz}, 1, 'last'))*c1*1000, P{yy, xx, zz}(find(P{yy, xx, zz}, 1):find(P{yy, xx, zz}, 1, 'last'))/Pmmax)
            %hold off
            % plot(P{yy, xx, zz}/ Pmmax)
            title(sprintf('X = %f, Y = %f, Z = %f', x(xx), y(yy), z(zz)));
            xlim([0, 4096])
            ylim([-2, 2])
            grid on
            grid minor
            pause(.1)
        end
    end
end