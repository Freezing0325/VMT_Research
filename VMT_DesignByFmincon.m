%% 初始配置

OriginStatus = 0;
% 目标序列
% GoalSequence = [1 0 0 1 0 1 1 1 0];
% GoalSequence = [1 1 0 0 0 1];
% RGB = true;
GoalSequence = [0 1 0];
RGB = false;
StepSum = size(GoalSequence, 2);
GoalSequence_Hat = [OriginStatus, GoalSequence(1: StepSum - 1)];
CalMethod = 2;
% 允许的最大归一化刚度
% Optimizer.MaxNormE = 9.320;
% Optimizer.MaxNormE = 4.691; %0.6粗负
% Optimizer.MaxNormE = 6.937; %0.6粗短
Optimizer.MaxNormE = 5;
% 允许的最小归一化刚度的差值，如果过小，在实际运行中，一侧的串联单元就不一定按照从小到大的顺序突跳
Optimizer.MinNormEDiff = 0.2;
Optimizer.MinNormEDiffFirstStage = 0.5;
% 允许的最小峰值点位置的差异，如果过小，在实际运行中，两侧的位移出现一定误差时就可能发生不同于设想的跳变，鲁棒性下降
Optimizer.MinDisDiff = 0.06;
% 对第一阶段峰值点位置差异的容许误差
Optimizer.MinDisDiffFirstStage = -0.01;
% 允许的最小的两个bit位之间位置的差异，如果过小，在实际运行中就可能出现两个切换的位置相互混淆的结果。
Optimizer.MinStepWall = 0.2;
% 允许的最大的结束时的位移差异，如果过大，那么在设计的序列切换结束后可能不稳定。
Optimizer.MaxOutDisDiff = 0;
% 允许的最大的突跳前力差异与输出单元突跳阈值之比，如果过于超过1，那么有可能在串联单元突跳前就使输出单元突跳至另一状态，或者在不需要突跳的时候发生突跳。
Optimizer.MaxFDiff = 0.95;

% 需要施加补偿的一侧，0：不需要，-1：左侧，1：右侧。
CompSide = GoalSequence - GoalSequence_Hat; 
if (~RGB)
    LeftComp = CompSide == -1;
    RightComp = CompSide == 1;
else
    LeftComp = CompSide ~= 1;
    RightComp = CompSide ~= -1;
    LeftComp(1) = 0;
    RightComp(1) = 0;
end

% CompSum(i)：前i个单元一共有几个补偿单元
CompSum = abs(CompSide);
for i = 2: StepSum
    CompSum(i) = CompSum(i - 1) + CompSum(i);
end
% 获得各个零势能点的位置
global Output_h a;
if (isempty(Output_h))
    VMT_Init();
end
OutputH = Output_h / a;
% OutputH = 0.0625;
U_0 = [0, (1: StepSum) - CompSum * OutputH * 2];  
BeginNormE = [1.5: 0.5: 1 + (StepSum-1) * 0.5, 1.5: 0.5: 1 + (StepSum-1) * 0.5];
[BestE, Bestg] = VMT_InverseDesign(GoalSequence, OriginStatus, CalMethod, BeginNormE, Optimizer);


%% 整理输出结果
CalMethod = 4;
LeftNormE = [1,BestE(1: StepSum - 1)];
RightNormE = [1,BestE(StepSum: 2 * (StepSum - 1))];
BestNormE = [LeftNormE, RightNormE];
FullNormE = BestNormE;
BestMaxNormE = max(BestNormE);

RealE = zeros(2, StepSum);
HeapPos = zeros(2, StepSum);
[RealE(1,:), HeapPos(1,:)] = VMT_CalHeapPos(LeftNormE, LeftComp, CalMethod);
[RealE(2,:), HeapPos(2,:)] = VMT_CalHeapPos(RightNormE, RightComp, CalMethod);

Judge_H = zeros(2, StepSum);
Delta_HeapPos = (OriginStatus - GoalSequence_Hat) * 2 * OutputH;
Judge_H(1,:) = double(HeapPos(1,:) + Delta_HeapPos);
Judge_H(2,:) = double(HeapPos(2,:) - Delta_HeapPos);
ChangeInfo = CompSide;

Real_H = zeros(2, StepSum);
Real_H(1,:) = double(Judge_H(1,:) + (ChangeInfo == -1) * 2 * OutputH);
Real_H(2,:) = double(Judge_H(2,:) + (ChangeInfo == 1) * 2 * OutputH);


fprintf('归一化刚度：\n');
for i = 1: 2 * StepSum
    fprintf('%.4f  ', BestNormE(i));
end
fprintf('\n');
fprintf('实际刚度：\n');
for i = 1: 2 * StepSum
    fprintf('%.4f  ', RealE(floor((i-1)/StepSum)+1, mod(i-1,StepSum)+1));
end
fprintf('\n');
if (~RGB)
    LeftComp = CompSide == -1;
    RightComp = CompSide == 1;
else
    LeftComp = CompSide ~= 1;
    RightComp = CompSide ~= -1;
    LeftComp(1) = 0;
    RightComp(1) = 0;
end
[PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE , LeftComp, RightComp, OriginStatus, CalMethod, []);
fprintf('预期序列：\n');
for i = 1: StepSum
    fprintf('%d  ', PredSequence(i));
end
fprintf('\n');
fprintf('最大力差异：%f\n', MaxForceDiff);

[MinDispDiff, MDDIndex] = min(abs(Judge_H(1,2:end)-Judge_H(2,2:end)));
fprintf('最小位移差异：%f，出现在第%d阶段\n', MinDispDiff, MDDIndex + 1);

global Normal_h
H_0 = Normal_h / a;
[Fm, ~] = VMT_SingleGetFm(1, H_0, CalMethod);
FinalDisDiff = (VMT_ConnectedGetU(RealE(1,:), H_0 - LeftComp * 2 * OutputH, BestMaxNormE * Fm, ones(1, StepSum), 2)...
                        - VMT_ConnectedGetU(RealE(2,:), H_0 - RightComp * 2 * OutputH, BestMaxNormE * Fm, ones(1, StepSum), 2)) * (1 - 2 * GoalSequence(StepSum));
fprintf('最终位移差异：%f\n', FinalDisDiff/OutputH); % 这个位移差异是考虑到最终状态时的结果，负值更稳定。

All_E = roundn(BestNormE, -4);
All_E_T = All_E.';

