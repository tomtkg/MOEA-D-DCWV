function p = Calp(points)
% Calculates the intermediate objective value p

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    [N,M] = size(points);
    Z = min(points,[],1);
    Zmax = max(points,[],1);
    Npoints = (points-repmat(Z,N,1))./repmat(Zmax-Z,N,1);

    normP   = sqrt(sum(Npoints.^2,2));
    CosineP = sum(Npoints./M,2).*sqrt(M)./normP;
    [~,I]   = min(normP.*sqrt(1-CosineP.^2));
    p       = normP(I)*CosineP(I) / sqrt(M);
end