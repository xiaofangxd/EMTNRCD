classdef NREG1 < PROBLEM
% <problem> <Sparse MOP>
% The feature selection problem
% dataNo --- 1 --- Number of dataset

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The datasets are taken from the UCI machine learning repository in
% http://archive.ics.uci.edu/ml/index.php
% No.   Name                              Samples Features Classes
% 1     MUSK1                               476     166       2
% 2     Semeion_handwritten_digit          1593     256      10
% 3     LSVT_voice_rehabilitation           126     310       2
% 4     ISOLET                             1557     617      26

    properties(Access = private)
        TrainIn;    % Input of training set
        TrainOut;   % Output of training set
    end
    properties(Access = public)
        labels;     % the labels of the node
    end
    methods
        %% Initialization
        function obj = NREG1()
            % Load data
            dataNo = obj.Global.ParameterSet(1);
            str = {'karate','polbooks','football','lesmis','dolphins','celegansneural'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_EGA_1.mat'),'DatasetA');
%             Data = Dataset.(str{dataNo});
            obj.TrainIn     = DatasetA.(str{dataNo});
            load(fullfile(fileparts(CallStack(1).file),'Dataset_EGY_1.mat'),'Datasety');
            obj.TrainOut    = Datasety.(str{dataNo});
            load(fullfile(fileparts(CallStack(1).file),[str{dataNo} '.mat']),'labels_real');
            obj.labels = labels_real;
            % Parameter setting
            obj.Global.M        = 2;
            obj.Global.D        = size(obj.TrainOut,2).^2;
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'binary';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec   = logical(PopDec);
            PopObj   = zeros(size(PopDec,1),obj.Global.M);
            for i = 1 : size(PopObj,1)
                TransferedPop = reshape(PopDec(i,:),sqrt(obj.Global.D),sqrt(obj.Global.D));
                % Clear self connecting
                TransferedPop(logical(eye(size(TransferedPop))))=0;
                for j = 1:sqrt(obj.Global.D)
                    temp = (obj.TrainIn(:,:,j)*TransferedPop(:,j)-obj.TrainOut(:,j));
                    PopObj(i,1) = sum(temp.^2)+PopObj(i,1);
                end
                PopObj(i,2) = sum(sum(TransferedPop));
%                 PopObj(i,3) = sum(sum(TransferedPop));
%                 %Community Detection
%                 if find(TransferedPop==1)
%                     modvec = cluster_jl(TransferedPop);
%                     %     modvec = cluster_jl_orient(varsafter);
%                     PopObj(i,2) = -modvec.MOD(end);
%                 else
%                     PopObj(i,2) = 0;
%                 end
            end
            PopObj(:,1) = PopObj(:,1)/(size(obj.TrainOut,1)*sqrt(obj.Global.D));
            
        end
    end
end