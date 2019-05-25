data = hdf5read('ref_pulse-40MHz.h5', 'ascan');
plot (data, 'ro')
hold on
plot (data, 'k')
hold off
grid on
grid minor