function CommunityNSGAII(Global)
% <algorithm> <W> 
%------------------------------- Copyright --------------------------------
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Last Update: 29.5.2022
%--------------------------------------------------------------------------

    %% Set the default parameters
    [psi,t1,q,delta] = Global.ParameterSet(3,1000,1,0.5);%%
    % The size of the population of weight-Individuals
    transformedProblemPopulationSize = 20;
    
    % There are different methods to select the xPrime solutions. 
    % Three methods can be chosen. The first one uses the largest Crowding
    % Distance values from the first non-dominated front. The second one
    % uses a tournament selection on the population based on
    % Pareto-dominance and Crowding Distance. The third option is 
    % introduced in publication (1), see above, based on
    % reference directions for the first m+1 chosen solutions and selects
    % random solutions afterwards. If method 3 is chosen, q is always at
    % least Global.M 
    methodToSelectxPrimeSolutions = 1;
    nmi = [];
    %% Generate random population
    Population = Global.Initialization();
    Global.NotTermination(Population,nmi);
    %% Optimise until end for uniformity. 
    remainingEvaluations    = Global.evaluation-delta*Global.evaluation;
    noOfParts               = floor(remainingEvaluations/1000);
%     
    for i = 1:noOfParts
        Population = fillPopulation(Population, Global);
        Population = Community_optimiseByNSGAII(Global, Population, 1000, false);
        Global.NotTermination(Population,nmi); 
    end
    xPrimeListfirst = Community_selectxPrimes(Population, q, methodToSelectxPrimeSolutions); 
%         [~, nmii] = Community_createGroups(xPrimeListfirst, Global.problem.labels);
    [Gfirst, nmii] = Community_createGroups_ECD1(xPrimeListfirst, Global.problem.labels);
%     [Gfirst, nmii] = Community_createGroups_DECS1(xPrimeListfirst, Global.problem.labels);
    nmi = [nmi nmii];
    %% Start the alternating optimisation 
    while Global.evaluated < Global.evaluation
        Population = Community_optimiseByNSGAII(Global, Population, Global.N, false);
        % Selection of xPrime solutions 
        xPrimeList = Community_selectxPrimes(Population, q, methodToSelectxPrimeSolutions); 
        WList   = [];
        % do for each xPrime
        for c = 1:size(xPrimeList,2)
            xPrime              = xPrimeList(c);
            % create variable groups
            % gamma----Number of nodes in network
            % G---The community index of each node
            % groups---[18,15,2]the number of nodes in community  
%             [G, nmii] = Community_createGroups(xPrime, Global.problem.labels); % changed
            [EP, G, nmii] = Community_createGroups_ECD2(xPrime, Gfirst, i, Global.problem.labels);
%                 [EP, G, nmii] = Community_createGroups_DECS2(xPrime, Gfirst, i, Global.problem.labels);
            nmi = [nmi nmii];
            Gfirst = G;
            groups = max(G);
            table=tabulate(G);
            [F,~]=max(table(:,2));
            I=find(table(:,2)==F);
            W2 = [];
            for indexgroups = 1:length(groups)+1
                % numberofoneC----Number of nodes in one community
                if indexgroups <= length(groups)
                    indexgroup = I(unidrnd(length(I)));
%                     indexgroup = unidrnd(groups);
                    eachgroups = find(G == indexgroup);
                else
                    eachgroups = [];
                end
                
                % a dummy object is needed to simulate the global class. Its
                % necessary to include this method into the Platemo framework.
                GlobalDummy = createGlobalDummy(xPrime, G, Global, transformedProblemPopulationSize, psi, eachgroups, indexgroups);
                
                % Create initial population for the transformed problem    
                DecomposePopulation = createInitialWeightPopulation(GlobalDummy);
                                
                % Optimise the transformed problem 
                DecomposePopulation = Community_optimiseByNSGAII(GlobalDummy, DecomposePopulation, t1-transformedProblemPopulationSize, true);

                % Extract the population 
                [W1, xPrime] = extractPopulation(DecomposePopulation, Global, Population, q, methodToSelectxPrimeSolutions, GlobalDummy, xPrime);
%                 W2 = INDIVIDUAL(W2);
                WList = [WList,W1];
            end
        end
        % Join populations. Duplicate solution (e.g. found in different
        % optimisation steps with different xPrimes) need to be removed. 
        Population          = eliminateDuplicates([Population,WList]);
        Population          = fillPopulation(Population, Global);
        
        % Environmental Selection
        [Population,~,~]    = Community_EnvironmentalSelection(Population,Global.N);
        if Global.evaluated >= Global.evaluation
            xPrimeListfirst = Community_selectxPrimes(Population, 1, 1); 
            [~, nmii] = Community_createGroups_nsga2_ECD1(xPrimeListfirst, Global.problem.labels);
            nmi(end) = nmii;
        end
        Global.NotTermination(Population,nmi);
    end
end

function GlobalDummy = createGlobalDummy(xPrime, G, Global, populationSize, psi, eachgroups, indexgroups)
    % Creates a dummy object. Needed to simulate the global class. Its
    % necessary to include this method into the Platemo
    % framework. 
    gamma = length(eachgroups);
    GlobalDummy = {};
    GlobalDummy.lower       = zeros(1,gamma);
    GlobalDummy.upper       = ones(1,gamma);
    GlobalDummy.N           = populationSize;
    GlobalDummy.xPrime      = xPrime;
    GlobalDummy.G           = G;
    GlobalDummy.psi         = psi;
    GlobalDummy.xPrimelower = Global.lower;
    GlobalDummy.xPrimeupper = Global.upper;
    GlobalDummy.isDummy     = true;
    GlobalDummy.Global      = Global;
    GlobalDummy.Index       = eachgroups;
    if isempty(eachgroups)
        tempD = sqrt(GlobalDummy.Global.D);
        x = reshape(GlobalDummy.xPrime.dec, tempD, tempD);
        for i = 1:length(unique(G))
            index = find(G == i);
            x(index,index) = inf;
        end
        xPrimeVar = reshape(x, 1, GlobalDummy.Global.D);
        [xIndex,yIndex] = find(xPrimeVar ~= inf);
        GlobalDummy.xIndex      = unique(xIndex);
        GlobalDummy.yIndex      = yIndex;
        GlobalDummy.D           = length(yIndex);
    else
        GlobalDummy.xIndex      = [];
        GlobalDummy.yIndex      = [];
        GlobalDummy.D           = gamma*gamma;
    end
    
    GlobalDummy.indexgroups = indexgroups;
end

function WeightPopulation = createInitialWeightPopulation(GlobalDummy) %% changed
    %creates an initial population for the transformed problem
    % N---the number of population
    WeightPopulation = [];
    decs = Community_WeightIndividual(GlobalDummy);
    for i = 1:GlobalDummy.N
        solution = C_WeightIndividual(decs(i,:),GlobalDummy);
        WeightPopulation = [WeightPopulation, solution];
    end
end

function Population = eliminateDuplicates(input)
    % Eliminates duplicates in the population
    [~,ia,~] = unique(input.objs,'rows');
    Population = input(ia);
end

function Population = fillPopulation(input, Global)
    % Fills the population with mutations in case its smaller than Global.N
    Population = input;
    theCurrentPopulationSize = size(input,2);
    if theCurrentPopulationSize < Global.N
        amountToFill    = Global.N-theCurrentPopulationSize;
        FrontNo         = NDSort(input.objs,inf);
        CrowdDis        = CrowdingDistance(input.objs,FrontNo);
        MatingPool      = TournamentSelection(2,amountToFill+1,FrontNo,-CrowdDis);
        Offspring       = GA(input(MatingPool));
        Population      = [Population,Offspring(1:amountToFill)];
    end
end

function [W1, W2] = extractPopulation(WeightPopulation, Global, Population, q, methodToSelectxPrimeSolutions, GlobalDummy, xPrime)
    % Extracts a population of individuals for the original problem based
    % on the optimised weights. 
    % First a selection of M+1 Weight-Individuals is selected and applied
    % to the whole population each. 
    % Second all Weight-Individuals are applied to the chosen xPrime
    % solution, since they are optimised for it. 
    weightIndividualList = Community_selectxPrimes(WeightPopulation, 1, 4); % need be changed
    PopDec1 = [];
    for wi = 1:size(WeightPopulation,2)
        PopDec1 = [PopDec1,WeightPopulation(wi).ind];
    end
    W1 = PopDec1;
    
    W2 = weightIndividualList.ind;
%     for wi = 1:size(WeightPopulation,2)
%         weightIndividual = WeightPopulation(wi);
%         if isempty(GlobalDummy.Index)
%             xIndex = GlobalDummy.xIndex;
%             yIndex = GlobalDummy.yIndex;
%             x = xP;
%             x(xIndex,yIndex) = weightIndividual.dec;         
%         else
%             weightVars = reshape(weightIndividual.dec, sqrt(GlobalDummy.D), sqrt(GlobalDummy.D));
%             x = reshape(xP, sqrt(Global.D), sqrt(Global.D));
%             x(V,V) = weightVars;
%             x = reshape(x, 1, Global.D);
%         end
%         PopDec2 = [PopDec2;x];
%     end
%     W3 = PopDec2;
end
