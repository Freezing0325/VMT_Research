%% 初始配置

OriginStatus = 1;
% 目标序列
% GoalSequence = [1 0 0 1 0 1 1 1 0];
GoalSequence = [1 0 1];
StepSum = size(GoalSequence, 2);

% 允许的最大归一化刚度
MaxNormE = 6;

% 允许的最小归一化刚度的差值，如果过小，在实际运行中，一侧的串联单元就不一定按照从小到大的顺序突跳
MinNormEDiff = 0.25;
% 允许的最小峰值点位置的差异，如果过小，在实际运行中，两侧的位移出现一定误差时就可能发生不同于设想的跳变，鲁棒性下降
MinDisDiff = 0.1;
% 允许的最小的两个bit位之间位置的差异，如果过小，在实际运行中就可能出现两个切换的位置相互混淆的结果。
MinStepWall = 0.2;
% 允许的最大的结束时的力差异，如果过大，那么在设计的序列切换结束后可能不稳定。
MaxOutDisDiff = 0.05;

% 迭代开始的归一化刚度配置
% BeginNormE = [1.5 2, 1.5 2];
% BeginNormE = [1.5     2       2.5     3         3.5     4       4.5     5, ...
%               1.5     2       2.5     3         3.5     4       4.5     5];
BeginNormE = [1.5: 0.5: 1 + (StepSum - 1) * 0.5, 1.5: 0.5: 1 + (StepSum - 1) * 0.5];

% 需要施加补偿的一侧，0：不需要，-1：左侧，1：右侧。
CompSide = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)]; 
LeftComp = CompSide == -1;
RightComp = CompSide == 1;

% CompSum(i)：前i个单元一共有几个补偿单元
CompSum = abs(CompSide);
for i = 2: StepSum
    CompSum(i) = CompSum(i - 1) + CompSum(i);
end
% 获得各个零势能点的位置
OutputH = 0.0625;
U_0 = [0, (1: StepSum) - CompSum * OutputH * 2];    


% 用来表示归一化刚度的符号矩阵
% TempNormE = sym('TempNormE_', [1, 2 * StepSum]);
% TempNormE(1) = 1;
% TempNormE(StepSum + 1) = 1;

% 用来表示所有待求归一化刚度的符号矩阵
NormE = sym('NormE_', [1, 2 * StepSum - 2]);
NormE_Mat = [1, NormE(1: StepSum - 1); 1, NormE(StepSum: 2 * (StepSum - 1))];
TempNormE = sym('TempNormE_', [1, 2 * StepSum]);

% 用符号表示的各个峰值位置，是归一化刚度的函数
[~, X_mL] = VMT_CalHeapPos_2(TempNormE(1: StepSum), LeftComp);
[~, X_mR] = VMT_CalHeapPos_2(TempNormE(StepSum + 1: 2 * StepSum), RightComp);

% A_SortE、B_SortE：要求归一化刚度是从小到大排列的，且相差至少为MinNormEDiff。
%
% A_SortE * NormE < B_SortE
%
% A_SortE_temp大致形状是：
% [-1    0   0   ...     0]
% [ 1   -1   0   ...     0]
% [ 0    1  -1   ...     0]
% [...  ... ...  ...    ..]
% [ 0    0   0   ...    -1]
% 这个矩阵组合成A_SortE：
% [A_SortE_temp     0       ]
% [     0       A_SortE_temp]
% B_SortE大致形状是：
% [-1 - M; -M; -M; ...; -M; -1 - M; -M; -M; ...; -M];

A_SortE_temp = [zeros(1, StepSum - 1); eye(StepSum - 2), zeros(StepSum - 2, 1)] - eye(StepSum - 1);
A_SortE = [A_SortE_temp, zeros(StepSum - 1); zeros(StepSum - 1), A_SortE_temp];
B_SortE = [-1; zeros(StepSum - 2, 1); -1; zeros(StepSum - 2, 1)] - MinNormEDiff * ones(2 * (StepSum - 1), 1);

global g_CallTimes con_CallTimes g_RunTime con_RunTime;
g_RunTime = 0;
con_RunTime = 0;
g_CallTimes = 0;
con_CallTimes = 0;
fprintf('开始迭代计算，过程可能要很久。\n')
AllRunTime = tic;
[BestE, Bestg] = fmincon(@(NormE)VMT_g_static([1, NormE(1: StepSum - 1), 1, NormE(StepSum: 2 * (StepSum - 1))], U_0, [X_mL; X_mR], GoalSequence, OriginStatus, MaxNormE), ...
                        BeginNormE, A_SortE, B_SortE, [], [], 1 * ones(1, 2 * (StepSum - 1)), MaxNormE * ones(1, 2 * (StepSum - 1)), ...
                        @(NormE)VMT_con_static([1, NormE(1: StepSum - 1), 1, NormE(StepSum: 2 * (StepSum - 1))], [X_mL; X_mR], GoalSequence, OriginStatus, U_0, MinDisDiff, MinStepWall, MaxNormE, MaxOutDisDiff));
toc(AllRunTime);

% QQ_Report('1603441246', 'Matlab算完了噢~');

%% 整理输出结果

[R_L, H_L] = VMT_CalHeapPos_2([1, BestE(1: StepSum - 1)], LeftComp);
[R_R, H_R] = VMT_CalHeapPos_2([1, BestE(StepSum: 2 * (StepSum - 1))], RightComp);

Delta_HeapPos = ~[OriginStatus, GoalSequence(1: StepSum - 1)] * 2 * OutputH;
Judge_H_L = double(H_L + Delta_HeapPos);
Judge_H_R = double(H_R - Delta_HeapPos);
Judge_H = [Judge_H_L; Judge_H_R];

Real_H_L = double(H_L + (~(GoalSequence & [OriginStatus, GoalSequence(1: StepSum - 1)])) * 2 * OutputH);
Real_H_R = double(H_R - (~(GoalSequence | [OriginStatus, GoalSequence(1: StepSum - 1)])) * 2 * OutputH);
Real_H = [Real_H_L; Real_H_R];

BestNormE = [1, BestE(1: StepSum - 1), 1, BestE(StepSum: 2 * (StepSum - 1))];
RealE = [R_L, R_R];
fprintf('归一化刚度：\n');
for i = 1: 2 * StepSum
    fprintf('%.4f  ', BestNormE(i));
end
fprintf('\n');
fprintf('实际刚度：\n');
for i = 1: 2 * StepSum
    fprintf('%.4f  ', RealE(i));
end
fprintf('\n');
LeftComp = CompSide == -1;
RightComp = CompSide == 1;
[PredSequence, MaxForceDiff] = VMT_GetSequence([1, BestE(1: StepSum - 1)], [1, BestE(StepSum: 2 * (StepSum - 1))] , LeftComp, RightComp, OriginStatus, 2, []);
fprintf('预期序列：\n');
for i = 1: StepSum
    fprintf('%d  ', PredSequence(i));
end
fprintf('\n');
fprintf('最大力差异：%f\n', MaxForceDiff);

All_E = roundn([1, BestE(1: StepSum - 1); 1, BestE(StepSum: 2 * (StepSum - 1))], -4);
All_E_T = All_E.';

