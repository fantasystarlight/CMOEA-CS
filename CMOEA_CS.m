classdef CMOEA_CS < ALGORITHM
% <multi/many> <real/binary/permutation><constrained/none>
% A comprehensive search based CMOEA
% V ---  0.1 --- The volume of potential list
% t ---  0.0001 --- The threshold of stage switching

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
            [V, t] = Algorithm.ParameterSet(0.1, 0.0001);
            [W, ~] = UniformPoint(floor(Problem.N/2), Problem.M);
            Stage = 1;
        %% Generate random population
            Population          = Problem.Initialization();
            FA = Population(unidrnd(Problem.N, [1, size(W, 2)]));
            DA = Population(unidrnd(Problem.N, [1, size(W, 2)]));
            FEA = Population;
            Pop1 = FA;
            Pop2 = DA; 
            PL = [];
            zmin = min(Population.objs,[],1) - 1e-6;    
        %% External archive initialization
            while Algorithm.NotTerminated(FEA)
        %% Optimization
                MatingPool_Pop1_1 = randperm(length(Pop1));
                MatingPool_Pop1_2 = randperm(length(Pop1));
                MatingPool_Pop2_1 = randperm(length(Pop2));
                MatingPool_Pop2_2 = randperm(length(Pop2));
                if rand() < 0.5
                    Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
                    Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
                else
                    Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
                    Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
                end
                Offspring = [Offspring1, Offspring2];
                if Stage == 1                    
                    zmin = min(zmin, min(Offspring.objs, [], 1) - 1e-6);  
                    [FA, FA_old] = ForwardExploreArchive(FA, Offspring, zmin, W);                   
                    DA = DiversityExploreArchive(DA, Offspring, zmin, W);
                    FEA = FeasibleExploitationArchive(FEA, Offspring, Problem.N);
                    Pop1 = FA;
                    Pop2 = DA;  
                    [PL, Stage] = SwitchStageStrategy(FA, FA_old, PL, V, t, Problem);
                else                    
                    zmin = min(FEA.objs, [], 1) - 1e-6;  
                    DA = DiversityExploreArchive(DA, Offspring, zmin, W);
                    FEA = FeasibleExploitationArchive(FEA, Offspring, Problem.N);
                    Pop1 = FEA(unidrnd(length(FEA),[1, size(W, 2)]));
                    Pop2 = DA;
                end
            end
        end
    end
end

