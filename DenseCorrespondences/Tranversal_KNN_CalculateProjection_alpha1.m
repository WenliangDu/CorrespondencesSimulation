function ProjectionSet = Tranversal_KNN_CalculateProjection_alpha1(Location1,TransMatrix)
%{
2015/12/16
Tranversal_KNN_CalculateProjection_alpha1
1. For calculating Projection Points by Location1 and TransMatrix
%}

Num = size(Location1,1);

ExtendedOne = ones(Num,1);
Location1Extended = [Location1 ExtendedOne];
ProjectionSetExtended =  TransMatrix * Location1Extended';


ProjectionSet = ( ProjectionSetExtended(1:2,:) )';