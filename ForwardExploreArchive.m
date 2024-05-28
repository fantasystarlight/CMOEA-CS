function [FA, FA_old] = ForwardExploreArchive(FA, Offspring, zmin, W)
% The Forward Exploration Archive of CMOEA-CS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu
    %% Dominance relation calculation
    FA_old = FA;
    [N, ~] = size(W);
    S = [FA, Offspring];
    FA = [];
    NonDominated = DominationCal(S, 0);
    S = S(NonDominated);
    Obj = S.objs;
    %% Enviornmental selection
    zmax = max(S.objs, [], 1); 
    Obj       = (Obj - zmin) ./ (zmax - zmin);
    Distance = sqrt(sum(Obj.^2, 2));
    Angle_S_to_W    = sin(acos(1 - pdist2(W,Obj,'cosine')));  
    for i = 1:N
        Fitness = Distance' .* Angle_S_to_W(i,:);
        [~, index] = min(Fitness);
        FA = [FA, S(index)];
    end   
end