function W = DCWVS(W, p)
% Distribution Control of Weight Vector Set

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [~,M] = size(W);
    A = W < 1.0 / M;
    W(A) = W(A) * p * M;
    W(~A) = 1.0 - (1.0 - W(~A)) * (1.0 - p) * M / (M - 1);
end