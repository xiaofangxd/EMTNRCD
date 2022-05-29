function [outIndexList,nmi] = Community_createGroups_DECS1(xPrime, labels)
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
    maxgen = 100;         % the maximum number of iterations
    pop_size = 100;       % the population size
    p_mutation = 0.20;    % the mutation rate
    p_migration = 0.50;   % the migration rate
    p_mu_mi = 0.50;       % the paramater to control the execution of mutation and migration
    PGLP_iter = 5;        % the number of iterations in PGLP 
    
    %community-based Grouping
    vars = xPrime.dec;
    vars(abs(vars)>0.05) = 1; vars(abs(vars)<0.05) = 0;
    % 1*MM->M*M
    numberOfNodes = sqrt(length(xPrime.dec));
    varsafter = reshape(vars,numberOfNodes,numberOfNodes);
    [dynMod, dynPop, outIndexList, dynTime] = DECS_1(varsafter, maxgen, pop_size, p_mutation, p_migration, p_mu_mi, PGLP_iter);
    nmi = NMI(labels, outIndexList);
%     nmii = NMI_t(labels, outIndexList);
%     modvec = cluster_jl(varsafter);
%     modvec = cluster_jl_orient(varsafter);
%     outIndexList = modvec.COM{end};
end