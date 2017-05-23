function [A_xy_array_Distorted,f] = TestWarping_DistortionFactors_alpha3(A_xy_array,x,y)
%{
2017/1/10
TestWarping_DistortionFactors_alpha1
1. Generate A_xy_array_Distorted by relative distortion.

2017/1/12
TestWarping_DistortionFactors_alpha2
1. Generate A_xy_array_Distorted by relative distortion for full-frame
(36mm*24mm), Radius=sqrt((36/2)^2 + (24/2)^2);

2017/01/15
TestWarping_DistortionFactors_alpha3
1. For centers = 0;
%}
%x = [0;5;15.25;23.15];y=[0;-0.5;-2.41;0]; % four points are good. %Distagon 2.8/21
%x = [0;10;15;20];y=[0;0.2;0.8;1.5]; % four points are good. %2.5 times ZEISS Otus 1.4/85
%x = [0;10;15;20];y=[0;0.4;1.6;3]; % four points are good. %2.5 times ZEISS Otus 1.4/85
%x = [0;10;17.5;23.15];y=[0;-1;-2;-1.8]; % four points are good. %2 times Distagon T* 2/25
%x = [0;10;17.5;23.15];y=[0;-2;-4;-3.6]; % four points are good. %2 times Distagon T* 2/25
[f,~] = TestMATLABFitting_alpha1(x,y);
%Radius=sqrt((36/2)^2 + (24/2)^2);
%Radius=sqrt((36/2)^2 + (24/2)^2);
%
%% DistortionFactors
u_o = [0 0];
R_origin = sqrt((A_xy_array(:,1)-u_o(1)).^2 + (A_xy_array(:,2)-u_o(2)).^2);
R_origin_Max = max(R_origin);
Miu = (R_origin./R_origin_Max).*23.15;
DistortionFactors = (Miu.^2).*f.a + (Miu.^4).*f.b;
DistortionFactors = DistortionFactors./100;

%% A_x_array_Distorted
A_x_array_Distorted = A_xy_array(:,1) - (A_xy_array(:,1) - u_o(1)).*DistortionFactors;
A_y_array_Distorted = A_xy_array(:,2) - (A_xy_array(:,2) - u_o(2)).*DistortionFactors;
A_xy_array_Distorted = [A_x_array_Distorted A_y_array_Distorted];



% % Plot Relative Distortion
% DistortionFactors_x = 0:0.05:23.15;
% DistortionFactors_y = (DistortionFactors_x.^2).*f.a + (DistortionFactors_x.^4).*f.b;
% figure,
% plot(DistortionFactors_x,DistortionFactors_y,'LineWidth',2');
% set(gca,'YTick',-2.5:0.5:0.5,'FontSize',18);
% set(gca,'XTick',0:5:25,'FontSize',18);
% set(gcf, 'position', [0 0 700 500]);
% 
% axis([0 25 -2.5 0.5]);
% grid on;
% xlabel('{Ideal Image Height(mm)}');
% ylabel('{Relative Distortion (%)}');