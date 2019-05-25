N=9000; %numero de pontos
py=120;
s2=load(['DADOS0H0.DAT']); % ou 78, ou 76, veja se tem muita diferença
s2_temp=reshape(s2,N,py);
vn = resample(s2_temp(:, 60), 320, 40);
vn = vn(find(vn == max(vn)):find(vn == max(vn)) + 1024)/max(abs(vn));
plot(vn)
k = numel(vn);
vn1 = vn(find(vn == max(vn)):end);
plot(vn)
hold on
plot(vn1)
vn1(end:end+10)=0;
plot(vn1)
vn1 = detrend(vn, 0);
plot(vn)


vhann = hann(find(vn1 == min(vn1)) + 800);
vhann(end:end+numel(vn1) - numel(vhann)) = 0;
plot(vn1)
hold on
plot(vhann)
plot(vhann.*vn1)
