function y = interface2(x, cr, df, dp, dpf)
y = x./sqrt(x.^2 + dp^2) - cr*(dpf - x)./sqrt((dpf - x).^2 + df^2);
end

% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.