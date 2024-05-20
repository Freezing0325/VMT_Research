function [PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, CompSide, OriginStatus, CalMethod, ActiveStatus)
%VMT_GetSequence       预测当前配置的VMT所对应的变化序列。
%
%   输入参数：
%   LeftNormE       左侧所有单元的归一化刚度
%   RightNormE      右侧所有单元的归一化刚度
%   CompSide        左右两侧单元的补偿情况，0：没有补偿，-1：左侧，1：右侧。
%   OriginStatus    初始状态
%   CalMethod       计算算法，1：求解非线性方程；2：线性模型：泰勒级数2阶拟合；3：非线性模型：泰勒级数3阶拟合。
%   ActiveStatus    两侧单元激活情况，第一行为左侧，第二行为右侧。不输入时，认为所有单元均激活。
%
%   输出：[PredSequence, MaxForceDiff]
%   PredSequence    预测的序列
%   MaxForceDiff    会影响序列变化的最大的力差异，用输出单元的突跳最大力得到比值

    OutputH = 0.0625;
    OutputL = 1.25;
    OutputE = 600;
    OutputFm = 2 * OutputE * 2/(3*sqrt(3)) * (OutputH/OutputL)^3 / sqrt((OutputH/OutputL)^2 + 1);
    Fm = 1/(6*sqrt(15)); 

    UnitSum = size(LeftNormE, 2);
    if (size(ActiveStatus, 1) == 0 || ~exist('ActiveStatus', 'var'))
        ActiveStatus = ones(2, UnitSum);
    end
    LeftActiveStatus = ActiveStatus(1, :);
    RightActiveStatus = ActiveStatus(2, :);
    LeftStepSum = sum(LeftActiveStatus ~= 0);
    RightStepSum = sum(RightActiveStatus ~= 0);
    if (LeftStepSum ~= RightStepSum)
        error('左右两侧激活数不同！\n');
    end
    StepSum = LeftStepSum;
    LeftComp = CompSide == -1;
    RightComp = CompSide == 1;
    [~, HeapPos_Left] = VMT_CalHeapPos(LeftNormE, LeftComp, CalMethod, LeftActiveStatus);
    [~, HeapPos_Right] = VMT_CalHeapPos(RightNormE, RightComp, CalMethod, RightActiveStatus);
    NowStatus = OriginStatus;
    PredSequence = zeros(1, StepSum);
    for i = 1: StepSum
        PredSequence(i) = HeapPos_Left(i) < HeapPos_Right(i) + 4 * OutputH * (NowStatus - 1);
        NowStatus = PredSequence(i);
    end

    CompInfo = [LeftComp; RightComp];
    ChangeInfo = PredSequence - [OriginStatus, PredSequence(1: StepSum - 1)];
    U_0 = zeros(2, StepSum);
    for i = 2: StepSum
        U_0(:, i) = U_0(:, i - 1) + 1 - CompInfo(:, i - 1) * OutputH * 4 + [-1; 1] * ChangeInfo(:, i - 1) * OutputH * 2;
    end
    
    MaxForceDiff = 0;
    % 在某阶段前的变形序列发生之后，这一步用来判断的峰值位置
    Delta_HeapPos = ~[OriginStatus, PredSequence(1: StepSum - 1)] * 2 * OutputH;
    Judge_HeapPos_L = HeapPos_Left + Delta_HeapPos;
    Judge_HeapPos_R = HeapPos_Right - Delta_HeapPos;
    % Real_HeapPos_L = HeapPos_Left + (~(PredSequence & [OriginStatus, PredSequence(1: StepSum - 1)])) * 2 * OutputH;
    % Real_HeapPos_R = HeapPos_Right - (~(PredSequence | [OriginStatus, PredSequence(1: StepSum - 1)])) * 2 * OutputH;
    for i = 1: StepSum
        if (ChangeInfo(i) == 0)
            continue;
        end
        if (PredSequence(i) == 1)
            % 防止出现预先突跳
            ThisDis = ((Judge_HeapPos_L(i) - U_0(2, i)) / (Judge_HeapPos_R(i) - U_0(2, i)) * RightNormE(i) - LeftNormE(i)) * Fm / OutputFm;
        else
            ThisDis = ((Judge_HeapPos_R(i) - U_0(1, i)) / (Judge_HeapPos_L(i) - U_0(1, i)) * LeftNormE(i) - RightNormE(i)) * Fm / OutputFm;
        end
        if (ThisDis > MaxForceDiff)
            MaxForceDiff = ThisDis;
        end
    end

end
