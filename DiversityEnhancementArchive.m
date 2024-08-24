function DA = DiversityEnhancementArchive(DA, Offspring, zmin, W)
% The Feasible Exploration Archive of CMOEA-CS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu
    %% Minimum angular distance calculation
    [N, ~] = size(W);
    Angle_W_to_W = acos(1 - pdist2(W,W,'cosine'));
    Angle_W_to_W(eye(N) == 1) = inf;
    h = mean(min(Angle_W_to_W));
    %% Dominance relation calculation
    S = [DA, Offspring];
    DA = [];
    S = DataCleaning(S, 1);
    Obj = S.objs;
    Con = S.cons;
    CV = sum(max(0,Con),2);
    %% Enviornmental selection
    if sum(CV <= 0) > 0
        zmax = max(Obj(CV <= 0, :), [], 1);
    else
        [~, index] = min(CV);
        zmax = Obj(index,:);
    end                    
    Obj       = (Obj - zmin) ./ (zmax - zmin);
    Angle_S_to_W    = sin(acos(1 - pdist2(W,Obj,'cosine')));
    for i = 1:N
        Angle = Angle_S_to_W(i,:);
        list = Angle <= h;
        if sum(list) == 0
            [~, index] = min(Angle_S_to_W(i,:));
            list(index) = true;
        end
        T = inf(length(S), 1);
        T(list) = CV(list);
        feasible = T <= 0;
        if sum(feasible) <= 0
            [~, index] = min(T);
        else
            T = inf(size(Angle));
            Fitness = Angle_S_to_W(i,:);
            T(feasible) = Fitness(feasible);
            [~, index] = min(T);
        end
        DA = [DA, S(index)];
    end        
end