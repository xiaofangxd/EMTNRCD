function [EP,outIndexList,nmi] = Community_createGroups_ECD2(xPrime, Gfirst, i, labels)
% Creates groups of the varibales. Three diffeent methods can be
% chosen. The first one uses linear groups, the second orders variables
% by absolute values, the third is a random grouping. For more
% information about these mechanisms see publication (2), see above. 
    
% ----------------------------------------------------------------------- 
%  WOF_createGroups.m 
%  Copyright (C) 2018 Heiner Zille
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
%  Author of this Code: 
%   Heiner Zille <heiner.zille@ovgu.de>
%
%  This file belongs to the following publications:
%
%  1) Heiner Zille and Sanaz Mostaghim
%     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
%     https://ieeexplore.ieee.org/document/8280974 
% 
%  2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
%     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
%     http://ieeexplore.ieee.org/document/7929324
%  
%  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
%     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
%     http://dl.acm.org/citation.cfm?id=2908979
%
%  Date of publication: 12.10.2018 
%  Last Update: 12.10.2018
% -----------------------------------------------------------------------
    % Parameter setting
    maxgen = 20;         % the maximum number of iterations
    pop_size = 100;       % the population size
    p_mutation = 0.20;    % the mutation rate
    p_migration = 0.50;   % the migration rate
    p_mu_mi = 0.50;       % the paramater to control the execution of mutation and migration
    Threshold = 0.80;     % R=1-Threshold is the parameter related to pupulation generation
    num_neighbor = 10;    % the neighbor size for each subproblem
    
    %community-based Grouping
    vars = xPrime.dec;
    vars(abs(vars)>0.05) = 1; vars(abs(vars)<0.05) = 0;
    % 1*MM->M*M
    numberOfNodes = sqrt(length(xPrime.dec));
    varsafter = reshape(vars,numberOfNodes,numberOfNodes);
    [dynMod, dynPop, outIndexList, EP, dynTime] = ECD_2(varsafter, maxgen, pop_size, p_mutation, p_migration, p_mu_mi, num_neighbor, Gfirst, Threshold, i);
    nmi = NMI(labels, outIndexList);
    if nmi>0.5
        www=1;
    end
    nmii = NMI_t(labels, outIndexList);
%     modvec = cluster_jl(varsafter);
%     modvec = cluster_jl_orient(varsafter);
%     outIndexList = modvec.COM{end};
%     groups = max(outIndexList);
end