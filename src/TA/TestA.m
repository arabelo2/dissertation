N=9000; %numero de pontos
py=120;
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H80.DAT']); % ou 78, ou 76, veja se tem muita diferença
s2_temp=reshape(s2,N,py);
%vn = resample(s2_temp(:, 60), 320, 40);
vn = s2_temp(:, 60);
%vn = resample(s2_temp(:, 60), 100, 100);
%vn = vn(find(vn == max(vn)) - 50:find(vn == max(vn)) + 125)/max(abs(vn));
plot(vn)
hold on
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H78.DAT']);
s2_temp=reshape(s2,N,py);
%vn = resample(s2_temp(:, 60), 320, 40);
vn = resample(s2_temp(:, 60), 100, 100);
%vn = vn(find(vn == max(vn)) - 50:find(vn == max(vn)) + 125)/max(abs(vn));
plot(vn)
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H76.DAT']);
s2_temp=reshape(s2,N,py);
%vn = resample(s2_temp(:, 60), 320, 40);
vn = resample(s2_temp(:, 60), 100, 100);
%vn = vn(find(vn == max(vn)) - 50:find(vn == max(vn)) + 125)/max(abs(vn));
plot(vn)

legend('DADOS0H80','DADOS0H78', 'DADOS0H76')