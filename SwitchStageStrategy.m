function [PL, Stage] = SwitchStageStrategy(FA, FA_old, PL, V, t, Problem)
% The stage switch strategy of CMOEA-CS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu
Volume = floor((Problem.maxFE / Problem.N) * V);
S = [FA, FA_old];
zmin = min(S.objs, [], 1) - 1e-6; 
zmax = max(S.objs, [], 1); 
FA_obj       = (FA.objs - zmin) ./ (zmax - zmin);
FA_old_obj     = (FA_old.objs - zmin) ./ (zmax - zmin);
Distance = sqrt(sum(FA_obj.^2, 2));
Distance_old = sqrt(sum(FA_old_obj.^2, 2));
MeanDistance = mean(Distance);
MeanDistance_old = mean(Distance_old);
PL = [PL, max(0, (MeanDistance_old - MeanDistance) / MeanDistance_old)];
Stage = 1;
if Problem.FE / Problem.maxFE > 0.8
    Stage = 2;
elseif length(PL) >= Volume
    if mean(PL) < t
        Stage = 2;
    else
        PL(1) = [];
    end
end
end