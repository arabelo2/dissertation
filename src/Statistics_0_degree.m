format longG;
xe = length(x_axis(x_axis <= 0));
xa = length(x(x <= 0));
% % Resample
picoapicointerpolado_x = resample(picoapico', 339, 100)'; size(picoapicointerpolado_x)
picoapicointerpolado_z = resample(picoapico, 3365, 1000); size(picoapicointerpolado_z)

[xds, yds] = datastats(picoapico(:,xe)/max(max(picoapico)), Ppp(xa, :)'/max(max(Ppp)))

% xds = 
% 
%        num: 201
%        max: 1
%        min: 0.321285140562249
%       mean: 0.559631561070152
%     median: 0.51004016064257
%      range: 0.678714859437751
%        std: 0.19754293474689
% 
% 
% yds = 
% 
%        num: 677
%        max: 1
%        min: 0.185349001376991
%       mean: 0.321917474131476
%     median: 0.236775673513049
%      range: 0.814650998623009
%        std: 0.201467625928193

[xds, yds] = datastats(picoapicointerpolado_z(:,xe)/max(max(picoapicointerpolado_z)), Ppp(xa, :)'/max(max(Ppp)))

% xds = 
% 
%        num: 677
%        max: 1
%        min: 0.0270568924038709
%       mean: 0.55758307594024
%     median: 0.509815081042639
%      range: 0.972943107596129
%        std: 0.198545000077201
% 
% 
% yds = 
% 
%        num: 677
%        max: 1
%        min: 0.185349001376991
%       mean: 0.321917474131476
%     median: 0.236775673513049
%      range: 0.814650998623009
%        std: 0.201467625928193

ze = length(z_axis(z_axis <= 40));
za = length(z(z <= 40/1000));
[xds, yds] = datastats(picoapico(ze,:)'/max(max(picoapico)), Ppp(:, za)/max(max(Ppp)))

% xds = 
% 
%        num: 120
%        max: 1
%        min: 0.0160642570281124
%       mean: 0.0669678714859437
%     median: 0.0441767068273092
%      range: 0.983935742971888
%        std: 0.130613735111218
% 
% 
% yds = 
% 
%        num: 407
%        max: 0.998756663537841
%        min: 0.0111230304229529
%       mean: 0.0416073103370638
%     median: 0.0274730398582235
%      range: 0.987633633114888
%        std: 0.0912640655599844

[xds, yds] = datastats(picoapicointerpolado_x(ze,:)'/max(max(picoapicointerpolado_x)), Ppp(:, za)/max(max(Ppp)))

xds = 
% 
%        num: 407
%        max: 1
%        min: 0.00303009212682293
%       mean: 0.0609820219370623
%     median: 0.03976138087115
%      range: 0.996969907873177
%        std: 0.11876273370223
% 
% 
% yds = 
% 
%        num: 407
%        max: 0.998756663537841
%        min: 0.0111230304229529
%       mean: 0.0416073103370638
%     median: 0.0274730398582235
%      range: 0.987633633114888
%        std: 0.0912640655599844

[xds, yds] = datastats(mag2db(picoapico(ze,:)'/max(max(picoapico))), mag2db(Ppp(:, za)/max(max(Ppp))))

% xds = 
% 
%        num: 120
%        max: 0
%        min: -35.8827871153555
%       mean: -27.1998157462497
%     median: -27.0961332387502
%      range: 35.8827871153555
%        std: 6.10335079035158
% 
% 
% yds = 
% 
%        num: 407
%        max: -0.0108062025590891
%        min: -39.0755374984829
%       mean: -30.9972138651506
%     median: -31.2218656753597
%      range: 39.0647312959238
%        std: 5.64258348072904

[xds, yds] = datastats(mag2db(picoapicointerpolado_x(ze,:)'/max(max(picoapicointerpolado_x))), mag2db(Ppp(:, za)/max(max(Ppp))))

% xds = 
% 
%        num: 407
%        max: 0
%        min: -50.3708833404318
%       mean: -28.054233152573
%     median: -28.0107708283984
%      range: 50.3708833404318
%        std: 6.19075222906115
% 
% 
% yds = 
% 
%        num: 407
%        max: -0.0108062025590891
%        min: -39.0755374984829
%       mean: -30.9972138651506
%     median: -31.2218656753597
%      range: 39.0647312959238
%        std: 5.64258348072904
