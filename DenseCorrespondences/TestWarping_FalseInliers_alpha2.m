function [SamplingNumbers,FalseInlierNumber,binoY,nomalS,SamplingL] = TestWarping_FalseInliers_alpha2(Ratio,ToltalL,A_xy_array)
%{
2017/01/10
TestWarping_FalseInliers_alpha1
1. Simulate FalseInliers as Additive white Gaussian noise.


2017/01/21
TestWarping_FalseInliers_alpha2
1. Remove: Ratio = Ratio/10;
%}
if Ratio == 0,
    SamplingNumbers = [];
    FalseInlierNumber = [];
    binoY = [];
    nomalS = [];
    SamplingL = [];
	return 
end
if isempty(Ratio),
    SamplingNumbers = [];
    FalseInlierNumber = [];
    binoY = [];
    nomalS = [];
    SamplingL = [];
	return 
end
%% normal distribution
SamplingL = round(ToltalL*Ratio);
SamplingNumbers = unidrnd(ToltalL,1,SamplingL);
SamplingNumbers = unique(SamplingNumbers);
SamplingNumbersL = length(SamplingNumbers);
% SamplingL = length(SamplingNumbers);
SamplingArray = 1:1:(ToltalL);
SamplingArray(SamplingNumbers) = [];
SamplingArrayL = length(SamplingArray);
while SamplingNumbersL ~= SamplingL,
    SamplingNumbersMore = unidrnd(SamplingArrayL,1,(SamplingL - SamplingNumbersL));
    SamplingNumbersMore = unique(SamplingNumbersMore);
    SamplingNumbers = [SamplingNumbers SamplingArray(SamplingNumbersMore)];
    
    %% Update
    SamplingArray(SamplingNumbersMore) = [];
    SamplingArrayL = length(SamplingArray);
    SamplingNumbersL = length(SamplingNumbers);
end
%SamplingNumbers = randi(ToltalL,1,SamplingL);

nomalS = 0:1:10;
binoY = binopdf(nomalS,10,1/2);

SelectNumber = round(binoY.*SamplingL);
if sum(SelectNumber) ~= SamplingL,
    SelectNumber(6) = SelectNumber(6) + (SamplingL - sum(SelectNumber));
end
SelectNumberL = length(SelectNumber);
%% Distance
LastNum = 0;
TraveralNum = zeros(SelectNumberL,2);
for i = 1:SelectNumberL, % TraveralNum
    if SelectNumber(i) ~= 0,
        TraveralNum(i,1) = LastNum + 1;
        TraveralNum(i,2) = TraveralNum(i,1) + SelectNumber(i) - 1;
        LastNum = TraveralNum(i,2);
    end
end
PercentageArray = zeros(SelectNumberL,2); % PercentageArray
PercentageTraverse = zeros(SelectNumberL,2);
InsertPercentage = 0;
for i = 1:SelectNumberL,
    %% PercentageArray
    PercentageArray(i,1) = InsertPercentage;
    if i == 1 || i == SelectNumberL,
        PercentageArray(i,2) = InsertPercentage + 5;
    else
        PercentageArray(i,2) = InsertPercentage + 10;
    end
    InsertPercentage = PercentageArray(i,2);
    
    %% PercentageTraverse
    if i == 1,
        PercentageTraverse(i,1) = 2;
    else
        PercentageTraverse(i,1) = round( (ToltalL/100)*PercentageArray(i,1) );
    end
    PercentageTraverse(i,2) = round( (ToltalL/100)*PercentageArray(i,2) );
    
end


SelectedRecord = true(1,ToltalL);
%SelectedRecord(SamplingNumbers) = false;

DistanceMatrix = sqrt((repmat(A_xy_array(:,1)',[SamplingL 1])-repmat(A_xy_array(SamplingNumbers,1),[1 ToltalL])).^2 + (repmat(A_xy_array(:,2)',[SamplingL 1])-repmat(A_xy_array(SamplingNumbers,2),[1 ToltalL])).^2);
FalseInlierNumber = zeros(1,SamplingL);
k = 1;
for i = 1:SelectNumberL,
    if SelectNumber(i) ~= 0,
        CurrentDistanceMatrix = DistanceMatrix(TraveralNum(i,1):TraveralNum(i,2),:);
        for j = 1:size(CurrentDistanceMatrix,1),
            [~,CurrentOrder] = sort(CurrentDistanceMatrix(j,:));
            CurrentTravereLeft = PercentageTraverse(i,1);
            CurrentTravereRight = PercentageTraverse(i,2);
            CurrentSelectOrderSub = CurrentOrder(CurrentTravereLeft:CurrentTravereRight);
            CurrentSelectOrder = CurrentSelectOrderSub(SelectedRecord(CurrentSelectOrderSub) == true);
            if ~isempty(CurrentSelectOrder),
                if length(CurrentSelectOrder) > 1,
                    CurrentSelectOrderSub2 = randi([1,length(CurrentSelectOrder)],1);
                else
                    CurrentSelectOrderSub2 = 1;
                end
            else
                CurrentSelectOrder = CurrentOrder(SelectedRecord(CurrentOrder) == true);
                CurrentSelectOrderSub2 = randi([1,length(CurrentSelectOrder)],1);
            end

            
            FalseInlierNumber(k) = CurrentSelectOrder(CurrentSelectOrderSub2);
            SelectedRecord(CurrentSelectOrder(CurrentSelectOrderSub2)) = false;
            k = k + 1;
        end
    end
end