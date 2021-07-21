function W = updateWeight(Objs,W)
% Update the weight vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    % Calculate the intermediate objective value p
    M = size(Objs,2);
    Objs = normalize(Objs,'range');
    normP   = sqrt(sum(Objs.^2,2));
    CosineP = sum(Objs./M,2).*sqrt(M)./normP;
    [~,I]   = min(normP.*sqrt(1-CosineP.^2));
    p       = normP(I)*CosineP(I) / sqrt(M);
    % Distribution control of weight vector set
    TF = W < 1.0 / M;
    W(TF) = W(TF) * p * M;
    W(~TF) = 1.0 - (1.0 - W(~TF)) * (1.0 - p) * M / (M - 1);
end