function Population = DataCleaning(Population, Add)
% The dominance based data cleaning of CMOEA-CS
% add = 0 --- calculate dominance relation without considering contraints
% add = 1 --- calculate dominance relation by considering contraints

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    PopObj = Population.objs;
    PopObj = roundn(PopObj, -10);
    [N, M] = size(PopObj);
    objs_dominated = false(1, N);
    PopCon = Population.cons;
    Cons = sum(max(0,PopCon),2);
    for i = 1: N-1
        err = PopObj(i,:) - PopObj;
        eq = zeros(N, 1);
        max_err = max(err, [],  2);
        min_err = min(err, [],  2);
        for j = i + 1: N
            for k = 1: M
                if err(j, k) ~= 0
                    break
                end
                if k == M
                    eq(j) = 1;
                end
            end
            if Add == 0
                if eq(j) == 1
                    objs_dominated(j) = true;
                elseif min_err(j) >= 0 
                    objs_dominated(i) = true;
                elseif max_err(j) <= 0
                    objs_dominated(j) = true;
                end
            else
                if eq(j) == 1
                    if Cons(i) <= Cons(j)
                        objs_dominated(j) = true;
                    else
                        objs_dominated(i) = true;
                    end
                elseif min_err(j) >= 0 
                    if Cons(j) <= 0 || Cons(j) <= Cons(i)
                        objs_dominated(i) = true;
                    end
                elseif max_err(j) <= 0
                    if Cons(i) <= 0 || Cons(i) <= Cons(j)
                        objs_dominated(j) = true;
                    end
                end
            end
        end
    end
    NonDominated = ~objs_dominated;
    Population = Population(NonDominated);
end