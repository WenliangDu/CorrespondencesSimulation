function [RandomConeCell,RealA_xy_array_3D,RealB_xy_array_3D,New_A_xy_array,New_B_xy_array,ConesPercentage] = TestWarping_RandomCone_alpha4(FocalLength,H,A_xy_array,RandomNum)
%{
2017/01/14
TestWarping_RandomCone_alpha1
1. Convert the script Test_Hight into a function.
2. Simulate the picture of the cones in the image plane.
3. Random cones as mountains or valley
Real means in realistic world

2017/01/21
TestWarping_RandomCone_alpha2
1. Calculate percentage of Cones

2017/01/29
TestWarping_RandomCone_alpha3
1. Use NewFactor

2017/05/23
TestWarping_RandomCone_alpha4
1. For Uploading
%}

EnLargeRatio = (H+FocalLength)/FocalLength;
TransMatrixScale = [EnLargeRatio 0 0;0 EnLargeRatio 0;0 0 1];


%% Generate RealA_xy_array (The positions of points in realistic.)
RealA_xy_array = Tranversal_KNN_CalculateProjection_alpha1(A_xy_array,TransMatrixScale);
RealCenter1 = [(max(RealA_xy_array(:,1))+min(RealA_xy_array(:,1)))/2 (max(RealA_xy_array(:,2))+min(RealA_xy_array(:,2)))/2];
UnitDistance = abs(RealA_xy_array(1,2) - RealA_xy_array(2,2)); % For overlapping Cones
%%
Overlapping_Array = RealA_xy_array( RealA_xy_array(:,1) >= RealCenter1(1,1) ,:);
Overlapping_Array_KDT = KDTreeSearcher(Overlapping_Array);
RealNum = 1:1:size(RealA_xy_array,1);
RealNumOfOverlapping = RealNum(RealA_xy_array(:,1) >= RealCenter1(1,1));
RecordOverlapping = true(1,size(Overlapping_Array,1));
%%
RealLl = max(RealA_xy_array(:,1));
RealLw = max(RealA_xy_array(:,2));
RealLwmin = min(RealA_xy_array(:,2));
%%
Quaterw = (max(RealA_xy_array(:,1))/2)*1;

RandomPosi = zeros(RandomNum,2);
MirroPosi = zeros(RandomNum,2); % For easy calculating
RandomConeCell = cell(RandomNum,4); 
% 1: ConeCenterNum(in the order of Overlapping_Array); 2: ConeCenterPosi; 
% 3: Hight of Cone; 4: Radius of Cone; 5: Points in Cone(in the order of
% Overlapping_Array); 6: Positions of Points in Cone; 7: Z of Points
%%
p=fittype('a*x+b','independent','x');
for i = 1:RandomNum, 
    %% Generate RandomPosi and MirroPosi
    RestOverlappingNum = RealNumOfOverlapping(RecordOverlapping);
    RandomConeCell{i,1} = RestOverlappingNum( randi([1,length(RestOverlappingNum)],1) );
    RandomConeCell{i,2} = RealA_xy_array(RandomConeCell{i,1},:);
    
    RandomPosi(i,:) = RandomConeCell{i,2};
    MirroPosi(i,1) = Quaterw - (RandomPosi(i,1) - Quaterw);
    MirroPosi(i,2) = RandomPosi(i,2);
    %%
    f1=fit([RealCenter1(1,1);RandomPosi(i,1)],[RealCenter1(1,2);RandomPosi(i,2)],p);
    f2=fit([RealCenter1(1,1);MirroPosi(i,1)],[RealCenter1(1,2);MirroPosi(i,2)],p);
    %% hMAX
    EndPosi1x = RealLl;
    EndPosi2x = RealLl;
    EndPosi1y = f1.a*EndPosi1x + f1.b;
    EndPosi2y = f2.a*EndPosi2x + f2.b;
    
    if f1.a > 0,
        if EndPosi1y > RealLw,
            EndPosi1y = RealLw;
            EndPosi1x = (EndPosi1y - f1.b)/f1.a;
        end
        if EndPosi2y > RealLw,
            EndPosi2y = RealLw;
            EndPosi2x = (EndPosi2y - f2.b)/f2.a;
        end
    else
        if EndPosi1y < RealLwmin,
            EndPosi1y = RealLwmin;
            EndPosi1x = (EndPosi1y - f1.b)/f1.a;
        end
        if EndPosi2y < RealLwmin,
            EndPosi2y = RealLwmin;
            EndPosi2x = (EndPosi2y - f2.b)/f2.a;
        end
    end
    
    Ditance1 = sqrt((RandomPosi(i,1) - EndPosi1x).^2 + (RandomPosi(i,2) - EndPosi1y).^2);
    Ditance2 = sqrt((MirroPosi(i,1) - EndPosi2x).^2 + (MirroPosi(i,2) - EndPosi2y).^2);
    
    if Ditance1 < Ditance2,
        Radius1 = sqrt((RealCenter1(1,1) - EndPosi1x).^2 + (RealCenter1(1,2) - EndPosi1y).^2);
        hMAX = (Ditance1/Radius1)*(FocalLength+H);
    else
        Radius2 = sqrt((RealCenter1(1,1) - EndPosi2x).^2 + (RealCenter1(1,2) - EndPosi2y).^2);
        hMAX = (Ditance2/Radius2)*(FocalLength+H);
    end
    
    %%
    Reasonalbe_h = hMAX/2; %
    alpha = (Reasonalbe_h)/(FocalLength+H);
    b1 = (alpha*sqrt((RealCenter1(1,1) - RandomPosi(i,1)).^2 + (RealCenter1(1,2) - RandomPosi(i,2)).^2))/(1-alpha);
    b2 = (alpha*sqrt((RealCenter1(1,1) - MirroPosi(i,1)).^2 + (RealCenter1(1,2) - MirroPosi(i,2)).^2))/(1-alpha);
    rMin = 2*( min(b1,b2) );

    PosiInCone = ( round(Overlapping_Array_KDT.X(:,1) - (RandomPosi(i,1)-rMin)) >= 0 & round(Overlapping_Array_KDT.X(:,1) - (RandomPosi(i,1)+rMin)) <= 0 ...
        & round(Overlapping_Array_KDT.X(:,2) - (RandomPosi(i,2)-rMin)) >= 0 & round(Overlapping_Array_KDT.X(:,2) - (RandomPosi(i,2)+rMin)) <= 0);
    CurrentRecordOverlapping = RecordOverlapping(PosiInCone);
    if ~all(CurrentRecordOverlapping), %% reduce rMin and Reasonalbe_h
        CoveredPosi = Overlapping_Array_KDT.X(PosiInCone,:);
        ProcessedPosi = CoveredPosi(~CurrentRecordOverlapping,:);

        %% 
        DistanceProcessed = sqrt((repmat(RandomPosi(i,1),[size(ProcessedPosi,1) 1]) - ProcessedPosi(:,1)).^2 + (repmat(RandomPosi(i,2),[size(ProcessedPosi,1) 1]) - ProcessedPosi(:,2)).^2);
        [minDis,minNum] = min(DistanceProcessed);
        if (ProcessedPosi(minNum,1) - RandomPosi(i,1)) == 0 || (ProcessedPosi(minNum,2) - RandomPosi(i,2)) == 0,
            New_rMin = minDis - UnitDistance;
        else
            New_rMin = (max( [abs(ProcessedPosi(minNum,1) - RandomPosi(i,1)) abs(ProcessedPosi(minNum,2) - RandomPosi(i,2))])) - UnitDistance;
        end
        if round(New_rMin) >= round(UnitDistance),
            PosiInCone = ( round(Overlapping_Array_KDT.X(:,1)-(RandomPosi(i,1)-New_rMin)) >= 0 & round(Overlapping_Array_KDT.X(:,1)-(RandomPosi(i,1)+New_rMin)) <= 0 ...
                & round(Overlapping_Array_KDT.X(:,2) - (RandomPosi(i,2)-New_rMin)) >= 0 & round(Overlapping_Array_KDT.X(:,2) - (RandomPosi(i,2)+New_rMin)) <= 0 );
            
            Reasonalbe_h = Reasonalbe_h*(New_rMin/rMin);
            rMin = New_rMin;
            CoveredNum = RealNumOfOverlapping(PosiInCone);
        else
            CoveredNum = RandomConeCell{i,1};
            Reasonalbe_h = Reasonalbe_h*((UnitDistance/2)/rMin);
            rMin = UnitDistance/2;
        end
    else
        CoveredNum = RealNumOfOverlapping(PosiInCone);
    end
    RandomConeCell{i,3} = Reasonalbe_h;
    RandomConeCell{i,4} = rMin;
    RandomConeCell{i,5} = CoveredNum;
    RandomConeCell{i,6} = RealA_xy_array(CoveredNum,:);
    RecordOverlapping(PosiInCone) = false;
    
    %% Calculate Z
    DistanceCone = sqrt((RandomConeCell{i,6}(:,1) - RandomConeCell{i,2}(1,1)).^2 + (RandomConeCell{i,6}(:,2) - RandomConeCell{i,2}(1,2)).^2);
    DistanceConeZ = ((repmat(rMin,size(DistanceCone)) - DistanceCone)./rMin).*Reasonalbe_h;
    DistanceConeZ(DistanceConeZ<0) = 0;
    
    RandomPositive = randi([0 1],1);
    if RandomPositive == 0,
        DistanceConeZ = DistanceConeZ.*(-1);
    end
    if isnan(DistanceConeZ),
        DistanceConeZ = 0;
    end
    RandomConeCell{i,7} = DistanceConeZ;
    
    %% Show Cones
%     figure(1),plot(RandomConeCell{i,6}(:,1),RandomConeCell{i,6}(:,2),'r*');
%     figure(2),plot3(RandomConeCell{i,6}(:,1),RandomConeCell{i,6}(:,2),RandomConeCell{i,7},'b*');hold on
    
end

%% RealA_xy_array_3D and RealB_xy_array_3D
RealA_xy_array_3D = RealA_xy_array;
RealA_xy_array_3D(:,3) = 0;
for j = 1:RandomNum,
    for k = 1:length(RandomConeCell{j,7})
        RealA_xy_array_3D(RandomConeCell{j,5}(k),3) = RandomConeCell{j,7}(k);
    end
end

RealB_xy_array_3D = RealA_xy_array;
RealB_xy_array_3D((RealB_xy_array_3D(:,1)<=0),3) = RealA_xy_array_3D((RealA_xy_array_3D(:,1)>=0),3);

%% New reference image plane: New_A_xy_array
New_A_xy_array = zeros(size(A_xy_array));
NewFactorUpper = (FocalLength + H);
NewFactor = repmat(NewFactorUpper,[size(RealA_xy_array_3D(:,3),1) 1])./(RealA_xy_array_3D(:,3).*(-1) + (FocalLength + H));

New_A_xy_array(:,1) = A_xy_array(:,1).*NewFactor;
New_A_xy_array(:,2) = A_xy_array(:,2).*NewFactor;
%% New template image plane: New_B_xy_array
New_B_xy_array = zeros(size(A_xy_array));
NewFactorB = repmat(NewFactorUpper,[size(RealB_xy_array_3D(:,3),1) 1])./(RealB_xy_array_3D(:,3).*(-1) + (FocalLength + H));
B_xy_array = A_xy_array;

New_B_xy_array(:,1) = B_xy_array(:,1).*NewFactorB;
New_B_xy_array(:,2) = B_xy_array(:,2).*NewFactorB;

%% Calculate Percentage
ConesPercentage = (length(find(RealA_xy_array_3D(RealA_xy_array_3D(:,3)~=0))) / (size(RealA_xy_array,1)/2))*100;
% figure,plot(B_xy_array(:,1),B_xy_array(:,2),'ro');hold on
% plot(New_B_xy_array(:,1),New_B_xy_array(:,2),'b*');
% for i = 1:length(New_B_xy_array),
%     if RealB_xy_array_3D(i,3) ~= 0,
%         line([B_xy_array(i,1) New_B_xy_array(i,1)],[B_xy_array(i,2) New_B_xy_array(i,2)]);
%     end
% end
% p=fittype('a*x.^2+b*x.^4','independent','x');
% f=fit(x,y,p);