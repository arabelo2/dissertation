N=9000; %numero de pontos
py=120;
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H0.DAT']); % ou 78, ou 76, veja se tem muita diferença
s2_temp=reshape(s2,N,py);
vn = resample(s2_temp(:, 60), 320, 100);
vn = vn(find(vn == max(vn))-215:find(vn == max(vn)) + 740)/max(abs(vn));
plot(vn)
grid on
grid minor
hold on
legend('DADOS0H0')
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H80.DAT']);
s2_temp=reshape(s2,N,py);
vn = resample(s2_temp(:, 60), 320, 100);
vn = vn(find(vn == max(vn))-215:find(vn == max(vn)) + 740)/max(abs(vn));
plot(vn)
grid on
grid minor
legend('DADOS0H80')
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H200.DAT']);
s2_temp=reshape(s2,N,py);
vn = resample(s2_temp(:, 60), 320, 100);
vn = vn(find(vn == max(vn))-215:find(vn == max(vn)) + 740)/max(abs(vn));
plot(vn)
grid on
grid minor
legend('DADOS0H200')

legend('DADOS0H0','DADOS0H80', 'DADOS0H200')