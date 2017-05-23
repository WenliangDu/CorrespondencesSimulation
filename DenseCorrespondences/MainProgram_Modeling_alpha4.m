
%{
2017/01/16
MainProgram_Modelling_alpha1
1. The main program of modelling camera

2017/2/2
MainProgram_Modeling_alpha3
1. Use ModelingCamera_alpha3

2017/5/23
MainProgram_Modeling_alpha4
1. For Uploading

Copyright (c) 2017 by Wenliang Du (du_wenliang@qq.com), 
Laboratory of Robotic Vision, School of Engineering Science, Macau University of Science and Technology
Last Modified October 2017/5/23
%}

%% Initialization the parameters of camera
FocalLength = 28; % 28 mm
CameraDistance = 40000;% 4 meters
%% Rigid parameters
% TranslateX = 50;
% TranslateY = 0;
% RotationAlpha = deg2rad(30);
TranslateX = 0;
TranslateY = 0;
RotationAlpha = 0;
%% Non-Rigid parameters
% ScaleT = 3;
% ShearPhi = deg2rad(30);
ScaleT = 1;
ShearPhi = 0;
%%
x = [0;10;15;20];y=[0;0.2;0.8;1.5]; % four points are good. %2.5 times ZEISS Otus 1.4/85
x1 = [0;10;17.5;23.15];y1=[0;-2;-4;-3.6]; % four points are good. %2 times Distagon T* 2/25
x0 = [0;5;15.25;23.15];y0=[0;-0.5;-2.41;0]; % four points are good. %Distagon 2.8/21
Distortion = cell(2,1);
Distortion{1,1} = [x,y];
Distortion{2,1} = [x1,y1];
Distortion{3,1} = [x0,y0];

%%
Ratio = 10;
Distortion_Sub = Distortion{1,1};
FeaturesNum = 1000;
RandomConesNum = 30;
%%
[MatchedLocation1,MatchedLocation2,TrueInliersIndex,FalseInliersIndex] = ModelingCamera_alpha4(Ratio,FeaturesNum,FocalLength,CameraDistance,RandomConesNum,TranslateX,TranslateY,RotationAlpha,ScaleT,ShearPhi,Distortion_Sub);







