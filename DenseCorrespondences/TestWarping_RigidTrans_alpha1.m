function B_xy_array = TestWarping_RigidTrans_alpha1(B_xy_array,TranslateX,TranslateY,RotationAlpha)
%{
2017/01/11
TestWarping_RigidTrans_alpha1
1. Rigid transforamtion.
%}
TransMatrixTranslate = [1 0 TranslateX;0 1 TranslateY;0 0 1];
TransMatrixRotation = [cos(RotationAlpha) sin(RotationAlpha) 0;-sin(RotationAlpha) cos(RotationAlpha) 0;0 0 1];
if TranslateX == 0 && TranslateY == 0 && RotationAlpha == 0,
    return
elseif TranslateX == 0 && TranslateY == 0,
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixRotation);
elseif RotationAlpha == 0,
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixTranslate);
else
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixRotation);
    B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixTranslate);
end





