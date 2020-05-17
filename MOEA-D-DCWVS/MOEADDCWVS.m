function MOEADDCWVS(Global)
% <algorithm> <M>
% Multiobjective evolutionary algorithm based on decomposition
% with Distribution Control of Weight Vector Set
% type --- 1 --- The type of aggregation function

%------------------------------- Reference --------------------------------
% Tomoaki Takagi, Keiki Takadama, and Hiroyuki Sato, A Distribution 
% Control of Weight Vector Set for Multi-objective Evolutionary Algorithms,
% Proc. of the Bio-inspired Information and Communication Technologies
% (BICT 2019), Lecture Notes of the Institute for Computer Sciences,
% Social Informatics and Telecommunications Engineering (LNICST),
% Vol. 289, Springer, Cham, pp. 70?80, 2019
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%----------------------------- Read me  Note ------------------------------
% If you want to use static approach, write main.m file  
% varargin = {'-algorithm',{@MOEADDCWVS, [A,p]}.
% If you want to use dynamic approach, write main.m file 
% varargin = {'-algorithm',{@MOEADDCWVS, A}.
% A specifies a scalarization function.
% p specifies the intermediate objective value in the range [0,1].
%--------------------------------------------------------------------------

% This algorithm is written by Tomoaki Takagi

    %% Parameter setting
    type = Global.ParameterSet(1);
    
    %% If there are no parameters, use the dynamic approach
    if size(type) == 1
         type = [type,0];
    end
    
    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);
    
    if type(2) == 0
        OriginalW = W;
    else
        % Static approach
        W = DCWVS(W,type(2));
    end    
    
    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        
        if type(2) == 0
            % Dynamic approach
            p = Calp(Population.objs);
            W = DCWVS(OriginalW,p);
        end
        
        % For each solution
        for i = 1 : Global.N      
            % Choose the parents
            P = B(i,randperm(size(B,2)));

            % Generate an offspring
            Offspring = GAhalf(Population(P(1:2)));

            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Update the neighbours
            switch type(1)
                case 1
                    % PBI approach
                    normW   = sqrt(sum(W(P,:).^2,2));
                    normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                    normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                    CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                    CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                    g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                case 2
                    % Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                case 3
                    % Tchebycheff approach with normalization
                    Zmax  = max(Population.objs,[],1);
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
                case 4
                    % Modified Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
                case 5
                    % Modified Tchebycheff approach with normalization
                    Zmax  = max(Population.objs,[],1);
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1)./W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1)./W(P,:),[],2);
            end
            Population(P(g_old>=g_new)) = Offspring;
        end
    end
end