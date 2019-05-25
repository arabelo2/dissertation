THETA = 20;
k = 2 * pi * f / c1;
ANGLE = -45:.001:45;
upsilon = 1/2 * k * B2x * (sin(deg2rad(ANGLE)) - sin(deg2rad(THETA))); 

figure()
H = abs(sin(upsilon)./upsilon); 
plot(ANGLE, 10*log10(H));
xlim([-45, 45]);
axis square
grid on
grid minor
xlabel('$${\theta (^{\circ})}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('10*log10(H)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel('Directivity $${(\theta)}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
legend(['{\theta_0 = }' num2str(THETA), '^{\circ}']);
set(gca,'FontSize',15);

figure()
H = abs(sin(upsilon)./upsilon); 
plot(ANGLE, H);
xlim([-45, 45]);
axis square
grid on
grid minor
xlabel('$${\theta (^{\circ})}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
ylabel(' Beam pattern $${(\theta)}$$ ', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
legend(['{\theta_0 = }' num2str(THETA), '^{\circ}']);
set(gca,'FontSize',15);