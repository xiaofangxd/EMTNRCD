function [Population,Dec,Mask] = Community_optimiseBySparseEA(GlobalDummy, inputPopulation, evaluations, isDummy, Dec, Mask, Fitness, REAL)
% <algorithm> <N>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = inputPopulation;
    [Population,Dec,Mask,FrontNo,CrowdDis] = Community_EnvironmentalSelection(Population,Dec,Mask,GlobalDummy.N);
	maximum = currentEvaluations(GlobalDummy, isDummy) + evaluations;
    %% Optimization
    while currentEvaluations(GlobalDummy, isDummy) < maximum
        MatingPool = TournamentSelection(2,GlobalDummy.N,FrontNo,-CrowdDis);
        [OffDec,OffMask] = Community_Operator(Dec(MatingPool,:),Mask(MatingPool,:),Fitness,REAL);
        if isDummy == true
            L = size(OffDec,1);
            Offspring = [];
            tep = OffDec.*OffMask;
            for i = 1:L
                Offspring = [Offspring, C_WeightIndividual(tep(i,:),GlobalDummy)];
            end 
        else
            Offspring = INDIVIDUAL(OffDec.*OffMask);
        end
        [Population,Dec,Mask,FrontNo,CrowdDis] = Community_EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],GlobalDummy.N);
    end
end

function e = currentEvaluations(GlobalDummy, isDummy)
    if isDummy == true  
        e = GlobalDummy.Global.evaluated;
    else
        e = GlobalDummy.evaluated;
    end
end