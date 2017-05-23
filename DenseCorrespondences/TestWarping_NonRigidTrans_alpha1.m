function B_xy_array = TestWarping_NonRigidTrans_alpha1(B_xy_array,ScaleT,ShearPhi)
%{
2017/01/11
TestWarping_NonRigidTrans_alpha1
1. NonRigid transforamtion.
%}
TransMatrixScale = [ScaleT 0 0;0 ScaleT 0;0 0 1];
TransMatrixShear = [1 tan(ShearPhi) 0;0 1 0;0 0 1];

if ScaleT == 1 && ShearPhi == 0,
    return
elseif ScaleT == 1,
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixShear);
elseif ShearPhi == 0,
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixScale);
else
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixScale);
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixShear);
end









