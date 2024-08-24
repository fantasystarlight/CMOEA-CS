classdef CMOEA_CS < ALGORITHM
% <multi/many> <real/binary/permutation><constrained/none>
% A comprehensive search based CMOEA

%------------------------------- Reference --------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    methods
        function main(Algorithm,Problem)
        %% Parameter setting
            NS = floor(Problem.N/3);
            [W, ~] = UniformPoint(NS, Problem.M);
        %% Generate random population
            Population          = Problem.Initialization();
            FA = [];
            DA = [];
            FEA = Population;
            Offspring = Population;
            zmin = min(Population.objs,[],1) - 1e-6;    
        %% External archive initialization
            while Algorithm.NotTerminated(FEA)
        %% Optimization                
                zmin = min(zmin, min(Offspring.objs, [], 1) - 1e-6);
                FA = ForwardExplorationArchive(FA, Offspring, zmin, W);                   
                DA = DiversityEnhancementArchive(DA, Offspring, zmin, W);
                FEA = FeasibleExploitationArchive(FEA, Offspring, Problem.N);
                Pop1 = FA;
                Pop2 = DA;  
                Pop3 = FEA(unidrnd(length(FEA), [1, NS]));
                MatingPool_Pop1_1 = randperm(length(Pop1));
                MatingPool_Pop1_2 = randperm(length(Pop1));
                MatingPool_Pop2_1 = randperm(length(Pop2));
                MatingPool_Pop2_2 = randperm(length(Pop2));
                MatingPool_Pop3_1 = randperm(length(Pop3));
                MatingPool_Pop3_2 = randperm(length(Pop3));
                if rand() < 0.5
                    Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
                    Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
                    Offspring3 = OperatorDE(Problem, Pop3, Pop3(MatingPool_Pop3_1), Pop3(MatingPool_Pop3_2),{1,0.5,1,1});
                else
                    Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
                    Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
                    Offspring3 = OperatorGA(Problem, Pop3(MatingPool_Pop3_1), {1,20,1,1});
                end
                Offspring = [Offspring1, Offspring2, Offspring3];                                                
            end
        end
    end
end

