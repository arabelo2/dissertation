xe = length(x_axis(x_axis <= 0));
xa = length(x(x <= 0));
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