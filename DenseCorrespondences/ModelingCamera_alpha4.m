function [Inliers_Reference,Inliers_Template,TrueInliersIndex,FalseInliersIndex,ConesPercentage] = ModelingCamera_alpha4(Ratio,FeaturesNum,FocalLength,CameraDistance,RandomConesNum,TranslateX,TranslateY,RotationAlpha,ScaleT,ShearPhi,Distortion)
%{
2017/1/10
TestWarping_alpha2
1. Simulate the camera distortion as a function.
BestRealRaito = [0.26 0.52 0.78 1.06 1.33 1.65 1.94 2.25 2.58 2.94];
IdealRatio = [5 10 15 20 25 30 35 40 45 50];

2017/01/14
TestWarping_alpha3
1. Add TestWarping_RandomCone_alpha1

2017/01/16
TestWarping_alpha4
1. Integrate TestWarping_RandomCone_alpha1 for random cones

2017/01/16
ModellingCamera_alpha1
1. Arrangement

2017/01/21
ModelingCamera_alpha2
1. Calculate percentage of Cones

2017/02/02
ModelingCamera_alpha3
1. Use TestWarping_RandomCone_alpha3

2017/05/23
ModelingCamera_alpha4
1. For uploading
%}

%% IdealRatio to Processed Ratio
% BestRealRaito = [0 0.26 0.52 0.78 1.06 1.33 1.65 1.94 2.25 2.58 2.94]; %line1 = 0:0.2:9;
IdealRatio = [0 10 20 30 40 50 60 70 80 90]; % (FalseLength/(TrueLength + FalseLength))*100;
BestRealRaito = [0 0.092 0.172 0.245 0.315 0.383 0.453 0.527 0.611 0.719]; %FeaturesNum = 5000;
ProcessedRatio = BestRealRaito(IdealRatio == Ratio);
warning('off');
%% Initial A_xy_array with Full-Frame Sensor (36mm*24mm /35mm)
FeaturesNum = FeaturesNum*2;
Lx = round(sqrt(FeaturesNum*36/24));
Ly = round(sqrt(FeaturesNum*24/36));
ToltalL = Lx*Ly;
if mod(ToltalL,2) ~= 0,
    Lx = Lx + 1;
    ToltalL = Lx*Ly;
end

linex = 1:1:Lx;
liney = 1:1:Ly;
A_x = repmat(linex,[Ly 1]);
A_y = repmat(flipud(liney'),[1 Lx]);

A_x_array = reshape(A_x,[ToltalL 1]);
A_y_array = reshape(A_y,[ToltalL 1]);
A_xy_array = [A_x_array A_y_array]; % Features in reference image


%% Move A_xy_array's center to zero
Center1 = [(max(A_xy_array(:,1))+min(A_xy_array(:,1)))/2 (max(A_xy_array(:,2))+min(A_xy_array(:,2)))/2];
TransMatrixTranslateToZero = [1 0 (0 - Center1(:,1));0 1 (0 - Center1(:,2));0 0 1];
A_xy_array = Tranversal_KNN_CalculateProjection_alpha1(A_xy_array,TransMatrixTranslateToZero);

[A_xy_array_Distorted,f] = TestWarping_DistortionFactors_alpha3(A_xy_array,Distortion(:,1),Distortion(:,2));
%% Random Cones
[~,~,~,New_A_xy_array,New_B_xy_array,ConesPercentage] = TestWarping_RandomCone_alpha4(FocalLength,CameraDistance,A_xy_array,RandomConesNum);

%% warping
[New_A_xy_array_Distorted,~] = TestWarping_DistortionFactors_alpha3(New_A_xy_array,Distortion(:,1),Distortion(:,2));
[New_B_xy_array_Distorted,~] = TestWarping_DistortionFactors_alpha3(New_B_xy_array,Distortion(:,1),Distortion(:,2));

%% Translation
New_A_xy_array_Distorted_RigidTransed = TestWarping_RigidTrans_alpha1(New_A_xy_array_Distorted,TranslateX,TranslateY,RotationAlpha);
New_A_xy_array_Distorted_NonRigidTransed = TestWarping_NonRigidTrans_alpha1(New_A_xy_array_Distorted_RigidTransed,ScaleT,ShearPhi);

%% Simulate FalseInliers as Additive white Gaussian noise
Overlapping_A_Index = (New_A_xy_array(:,1) >= 0);
Overlapping_B_Index = (New_B_xy_array(:,1) <= 0);
Overlapping_New_A_xy_array = New_A_xy_array(Overlapping_A_Index,:);
Overlapping_A_xy_array = A_xy_array(New_A_xy_array(:,1) >= 0,:);
Overlapping_ToltalL = size(Overlapping_New_A_xy_array,1);

[SamplingNumbers,FalseInlierNumber,binoY,nomalS,SamplingL] = TestWarping_FalseInliers_alpha2(ProcessedRatio,Overlapping_ToltalL,Overlapping_A_xy_array);

RecordTrueInliers = true(1,Overlapping_ToltalL);
RecordTrueInliers(SamplingNumbers) = false;
RecordTrueInliers(FalseInlierNumber) = false;
RecordOverlappingNum = 1:1:Overlapping_ToltalL;
TrueInliersNumbers = RecordOverlappingNum(RecordTrueInliers);

%% Overlapping area
Overlapping_New_B_xy_array_Distorted = New_B_xy_array_Distorted(Overlapping_B_Index,:);
Overlapping_New_A_xy_array_Distorted_NonRigidTransed = New_A_xy_array_Distorted_NonRigidTransed(Overlapping_A_Index,:);

%% Outputs
Inliers_Reference = [Overlapping_New_A_xy_array_Distorted_NonRigidTransed(TrueInliersNumbers,:) ; Overlapping_New_A_xy_array_Distorted_NonRigidTransed(SamplingNumbers,:)];
Inliers_Template = [Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers,:) ; Overlapping_New_B_xy_array_Distorted(FalseInlierNumber,:)];
TrueInliersIndex = 1:1:length(TrueInliersNumbers);
FalseInliersIndex = (length(TrueInliersNumbers) + 1):1:size(Inliers_Reference,1);

%% Figures
%{
figure,
plot(A_xy_array_Distorted(:,1),A_xy_array_Distorted(:,2),'r*'); hold on
plot(A_xy_array(:,1),A_xy_array(:,2),'g*');


%
Overlapping_New_A_xy_array_Distorted = New_A_xy_array_Distorted(Overlapping_A_Index,:);
TempTransMatrix = Transformation_alpha1(Overlapping_New_A_xy_array_Distorted_NonRigidTransed,Overlapping_New_B_xy_array_Distorted);
Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected = Tranversal_KNN_CalculateProjection_alpha1(Overlapping_New_A_xy_array_Distorted_NonRigidTransed,TempTransMatrix);
%% Show Distorted
CenterA_Distorted = [mean(Overlapping_New_A_xy_array_Distorted(:,1)) mean(Overlapping_New_A_xy_array_Distorted(:,2))];
CenterB_Distorted = [mean(Overlapping_New_B_xy_array_Distorted(:,1)) mean(Overlapping_New_B_xy_array_Distorted(:,2))];
TransMatrixlateAtoB_Overlapping_Distorted = [1 0 (CenterB_Distorted(:,1) - CenterA_Distorted(:,1));0 1 (CenterB_Distorted(:,2) - CenterA_Distorted(:,2));0 0 1];
Translated_Overlapping_New_A_xy_array_Distorted = Tranversal_KNN_CalculateProjection_alpha1(Overlapping_New_A_xy_array_Distorted,TransMatrixlateAtoB_Overlapping_Distorted);
figure,
plot(Translated_Overlapping_New_A_xy_array_Distorted(:,1),Translated_Overlapping_New_A_xy_array_Distorted(:,2),'ro');hold on
plot(Overlapping_New_B_xy_array_Distorted(:,1),Overlapping_New_B_xy_array_Distorted(:,2),'b*');
plot(Translated_Overlapping_New_A_xy_array_Distorted(TrueInliersNumbers,1),Translated_Overlapping_New_A_xy_array_Distorted(TrueInliersNumbers,2),'mo');
plot(Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers,1),Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers,2),'c*');
if ~isempty(SamplingNumbers),
    plot(Translated_Overlapping_New_A_xy_array_Distorted(SamplingNumbers,1),Translated_Overlapping_New_A_xy_array_Distorted(SamplingNumbers,2),'ro');
    plot(Overlapping_New_B_xy_array_Distorted(FalseInlierNumber,1),Overlapping_New_B_xy_array_Distorted(FalseInlierNumber,2),'b*');
    %% Draw lines
    for i = 1:length(TrueInliersNumbers),
        line([Translated_Overlapping_New_A_xy_array_Distorted(TrueInliersNumbers(i),1) Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers(i),1)],[Translated_Overlapping_New_A_xy_array_Distorted(TrueInliersNumbers(i),2) Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers(i),2)],'Color','g','LineWidth',1)
    end
    for i = 1:length(SamplingNumbers),
        line([Translated_Overlapping_New_A_xy_array_Distorted(SamplingNumbers(i),1) Overlapping_New_B_xy_array_Distorted(FalseInlierNumber(i),1)],[Translated_Overlapping_New_A_xy_array_Distorted(SamplingNumbers(i),2) Overlapping_New_B_xy_array_Distorted(FalseInlierNumber(i),2)],'Color','r','LineWidth',1)
    end

end
%% Show Transfermation
figure,
plot(New_A_xy_array_Distorted(:,1),New_A_xy_array_Distorted(:,2),'ro'); hold on
plot(New_A_xy_array_Distorted_NonRigidTransed(:,1),New_A_xy_array_Distorted_NonRigidTransed(:,2),'b*');

figure,
plot(Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(:,1),Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(:,2),'ro'); hold on
plot(Overlapping_New_B_xy_array_Distorted(:,1),Overlapping_New_B_xy_array_Distorted(:,2),'b*');
plot(Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(TrueInliersNumbers,1),Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(TrueInliersNumbers,2),'mo');
plot(Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers,1),Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers,2),'c*');
if ~isempty(SamplingNumbers),
    plot(Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(SamplingNumbers,1),Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(SamplingNumbers,2),'ro');
    plot(Overlapping_New_B_xy_array_Distorted(FalseInlierNumber,1),Overlapping_New_B_xy_array_Distorted(FalseInlierNumber,2),'b*');
    %% Draw lines
    for i = 1:length(TrueInliersNumbers),
        line([Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(TrueInliersNumbers(i),1) Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers(i),1)],[Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(TrueInliersNumbers(i),2) Overlapping_New_B_xy_array_Distorted(TrueInliersNumbers(i),2)],'Color','g','LineWidth',1)
    end
    for i = 1:length(SamplingNumbers),
        line([Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(SamplingNumbers(i),1) Overlapping_New_B_xy_array_Distorted(FalseInlierNumber(i),1)],[Overlapping_New_A_xy_array_Distorted_NonRigidTransed_Projected(SamplingNumbers(i),2) Overlapping_New_B_xy_array_Distorted(FalseInlierNumber(i),2)],'Color','r','LineWidth',1)
    end
end

%% Show Outputs
figure,
plot(Inliers_Reference(TrueInliersIndex,1),Inliers_Reference(TrueInliersIndex,2),'mo');hold on
plot(Inliers_Template(TrueInliersIndex,1),Inliers_Template(TrueInliersIndex,2),'c*');

plot(Inliers_Reference(FalseInliersIndex,1),Inliers_Reference(FalseInliersIndex,2),'ro');
plot(Inliers_Template(FalseInliersIndex,1),Inliers_Template(FalseInliersIndex,2),'b*');
% Draw lines
for i = 1:length(TrueInliersIndex),
    line([Inliers_Reference(TrueInliersIndex(i),1) Inliers_Template(TrueInliersIndex(i),1)],[Inliers_Reference(TrueInliersIndex(i),2) Inliers_Template(TrueInliersIndex(i),2)],'Color','g','LineWidth',1)
end
for i = 1:length(FalseInliersIndex),
    line([Inliers_Reference(FalseInliersIndex(i),1) Inliers_Template(FalseInliersIndex(i),1)],[Inliers_Reference(FalseInliersIndex(i),2) Inliers_Template(FalseInliersIndex(i),2)],'Color','r','LineWidth',1)
end


%}
%% Show overlapping 
% B_xy_array = A_xy_array;
Overlapping_A_Indexx = find(A_xy_array(:,1) >= 0);
Overlapping_B_Indexx = find(A_xy_array(:,1) <= 0);
% 
TransMatrixx = [1 0 max(A_xy_array(:,1));0 1 0;0 0 1];
% B_xy_array = Tranversal_KNN_CalculateProjection_alpha1(B_xy_array,TransMatrixx);
% figure,
% plot(A_xy_array(:,1),A_xy_array(:,2),'r.','MarkerSize',8);hold on
% plot(B_xy_array(:,1),B_xy_array(:,2),'b*','MarkerSize',6);
% %axis([-10 20 -6.5 6.5]);
% set(gcf, 'position', [0 0 800 500]);
%% outliers
New_B_xy_array_Distorted1 = Tranversal_KNN_CalculateProjection_alpha1(New_B_xy_array_Distorted,TransMatrixx);
Overlapping_A_xy_array = New_A_xy_array_Distorted(Overlapping_A_Indexx,:);
Overlapping_B_xy_array = New_B_xy_array_Distorted1(Overlapping_B_Indexx,:);
figure,
plot(Overlapping_A_xy_array(:,1),Overlapping_A_xy_array(:,2),'ro','MarkerSize',6);hold on
plot(Overlapping_B_xy_array(:,1),Overlapping_B_xy_array(:,2),'b*','MarkerSize',6);
for i = 1:length(SamplingNumbers),
    line([Overlapping_A_xy_array(SamplingNumbers(i),1) Overlapping_B_xy_array(FalseInlierNumber(i),1)],[Overlapping_A_xy_array(SamplingNumbers(i),2) Overlapping_B_xy_array(FalseInlierNumber(i),2)],'Color','m','LineWidth',1)
end
for i = 1:length(TrueInliersNumbers),
    line([Overlapping_A_xy_array(TrueInliersNumbers(i),1) Overlapping_B_xy_array(TrueInliersNumbers(i),1)],[Overlapping_A_xy_array(TrueInliersNumbers(i),2) Overlapping_B_xy_array(TrueInliersNumbers(i),2)],'Color','g','LineWidth',2)
end
set(gcf, 'position', [0 0 800 500]);
title('Dense correspodences with random outliers');
%axis([-1 9 -7 7])

%% Show cones
% New_B_xy_array1 = Tranversal_KNN_CalculateProjection_alpha1(New_B_xy_array,TransMatrixx);
% figure,
% plot(New_A_xy_array(:,1),New_A_xy_array(:,2),'ro','MarkerSize',12);hold on
% plot(New_B_xy_array1(:,1),New_B_xy_array1(:,2),'bs','MarkerSize',12);
% axis([-10 20 -6.5 6.5]);
% set(gcf, 'position', [0 0 800 500]);

%% Show Distortion and cones
% New_B_xy_array_Distorted1 = Tranversal_KNN_CalculateProjection_alpha1(New_B_xy_array_Distorted,TransMatrixx);
% figure,
% plot(New_A_xy_array_Distorted(:,1),New_A_xy_array_Distorted(:,2),'ro','MarkerSize',6);hold on
% plot(New_B_xy_array_Distorted1(:,1),New_B_xy_array_Distorted1(:,2),'b*','MarkerSize',6);
% %axis([-10 20 -6.5 6.5]);
% set(gcf, 'position', [0 0 800 500]);
% title('feature points');
%% Show outliers
% figure,
% plot(New_A_xy_array_Distorted(:,1),New_A_xy_array_Distorted(:,2),'r.','MarkerSize',8);hold on
% plot(New_B_xy_array_Distorted1(:,1),New_B_xy_array_Distorted1(:,2),'b*','MarkerSize',6);
% 
% for i = 1:length(SamplingNumbers),
%     line([New_A_xy_array_Distorted(Overlapping_A_Indexx(SamplingNumbers(i)),1) New_B_xy_array_Distorted1(Overlapping_B_Indexx(FalseInlierNumber(i)),1)],[New_A_xy_array_Distorted(Overlapping_A_Indexx(SamplingNumbers(i)),2) New_B_xy_array_Distorted1(Overlapping_B_Indexx(FalseInlierNumber(i)),2)],'Color','m','LineWidth',0.5)
% end
% %TrueInliersNumbersIndex = find(TrueInliersNumbers);
% for i = 1:length(TrueInliersNumbers),
%     line([New_A_xy_array_Distorted(Overlapping_A_Indexx(TrueInliersNumbers(i)),1) New_B_xy_array_Distorted1(Overlapping_B_Indexx(TrueInliersNumbers(i)),1)],[New_A_xy_array_Distorted(Overlapping_A_Indexx(TrueInliersNumbers(i)),2) New_B_xy_array_Distorted1(Overlapping_B_Indexx(TrueInliersNumbers(i)),2)],'Color',[0,0.784,0],'LineWidth',2)
% end
% %axis([-10 20 -6.5 6.5]);
% set(gcf, 'position', [0 0 800 500]);


%% Plot Binomial Distribution
% if ~isempty(binoY),
%     figure,
%     binoYBar = bar(binoY,'r');
%     ch = get(binoYBar,'children');
%     %set(gca,'XTickLabel',{'0%~5%','5%~15%','15%~25%','25%~35%','35%~45%','45%~55%','55%~65%','65%~75%','75%~85%','85%~95%','95%~100%'})
%     set(gca,'XTickLabel',{'5%','10%','20%','30%','40%','50%','60%','70%','80%','90%','95%'},'FontSize',22);
%     set(gca,'YTick',0:0.05:0.26,'FontSize',26);
%     axis([0 12 0.0 0.26]) 
%     ylabel('Ratio of Sampling Number','FontSize',26);
%     xlabel('Outliers Ratio','FontSize',26);
%     grid on;
%     set(gcf, 'position', [0 0 1200 700]);
% %     for i = 1:length(binoY),
% %         text(nomalS(i)+0.7,binoY(i)+0.005,num2str(roundn(binoY(i),-4)));
% %     end
% end
%% Plot Relative Distortion
% DistortionFactors_x = 0:0.05:23.15;
% DistortionFactors_y = (DistortionFactors_x.^2).*f.a + (DistortionFactors_x.^4).*f.b;
% figure,
% plot(DistortionFactors_x,DistortionFactors_y,'LineWidth',2');
% set(gca,'YTick',-2.5:0.5:0.5,'FontSize',22);
% set(gca,'XTick',0:5:25,'FontSize',22);
% set(gcf, 'position', [0 0 700 500]);
% 
% axis([0 25 -2.5 0.5]);
% grid on;
% xlabel('{Ideal Image Height(mm)}','FontSize',22);
% ylabel('{Relative Distortion (%)}','FontSize',22);

