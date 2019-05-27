xe = length(x_axis(x_axis <= 40*sin(20*pi/180)));
xa = length(x(x <= 40*sin(20*pi/180)/1000));
[xds, yds] = datastats(picoapico(:,xe)/max(max(picoapico)), Ppp(xa, :)'/max(max(Ppp)))

% xds = 
% 
%        num: 201
%        max: 0.898305084745763
%        min: 0.0338983050847458
%       mean: 0.23273463192512
%     median: 0.169491525423729
%      range: 0.864406779661017
%        std: 0.194058517050416
% 
% 
% yds = 
% 
%        num: 677
%        max: 0.998819218283403
%        min: 0.0254282724099457
%       mean: 0.140569489603864
%     median: 0.0886135754927694
%      range: 0.973390945873457
%        std: 0.161128982110616


ze = length(z_axis(z_axis <= 40*cos(20*pi/180)));
za = length(z(z <= 40*cos(20*pi/180)/1000));
[xds, yds] = datastats(picoapico(ze,:)'/max(max(picoapico)), Ppp(:, za)/max(max(Ppp)))

% xds = 
% 
%        num: 120
%        max: 0.983050847457627
%        min: 0.0254237288135593
%       mean: 0.136228813559322
%     median: 0.0932203389830508
%      range: 0.957627118644068
%        std: 0.146842185671639
% 
% 
% yds = 
% 
%        num: 407
%        max: 0.998819218283403
%        min: 0.0106285235022538
%       mean: 0.0748434332054241
%     median: 0.0500638265216191
%      range: 0.988190694781149
%        std: 0.116144917674449