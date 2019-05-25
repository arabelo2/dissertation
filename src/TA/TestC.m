N=9000; %numero de pontos
py=120;
s2=load(['DADOS0H0.DAT']); % ou 78, ou 76, veja se tem muita diferença
s2_temp=reshape(s2,N,py);
vn = resample(s2_temp(:, 60), 320, 40);
k = numel(vn);
vn1 = vn(find(vn == max(vn))-10:end);
plot(vn)
hold on
plot(vn1)
vn1(end:end+k - numel(vn1))=0;
plot(vn1)
vn1 = detrend(vn, 0);
plot(vn)
