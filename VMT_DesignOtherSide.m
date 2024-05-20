function [BestRealE_Left, BestE, RealE_Right] = VMT_DesignOtherSide(RightNormE, GoalSequence, OriginStatus, BeginE)
%VMT_DesignOtherSide        根据已有的一侧设计与目标序列，设计另一侧的序列。
%
%   输入：RightNormE, GoalSequence, OriginStatus, [BeginE]
%   RightNormE：     另一侧的归一化刚度
%   GoalSequence：   目标序列
%   OriginStatus：   初始状态
%   [BeginE]：       迭代起始的归一化刚度，如果不输入则默认随机开始迭代。
%
%   输出：[BestRealE_Left, BestE, RealE_Right]
%   BestRealE_Left： 设计一侧的真实刚度
%   BestE：          设计一侧的归一化刚度
%   RealE_Right      另一侧的真实刚度
CalMethod = 2;

StepSum = size(GoalSequence, 2);    % 序列长度

OutputH = 0.0625;
OutputL = 1.25;
OutputE = 600;
OutputFm = 2 * OutputE * 2/(3*sqrt(3)) * (OutputH/OutputL)^3 / sqrt((OutputH/OutputL)^2 + 1);
Fm = 1/(6*sqrt(15)); 

IfChanged = [0, GoalSequence(1: end - 1) ~= OriginStatus];  % 目前的状态是否不同于初始状态
CompSide = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)]; % 需要施加补偿的一侧，0：不需要，-1：左侧，1：右侧。
LeftComp = CompSide == -1;
RightComp = CompSide == 1;
[RealE_Right, HeapPos_Right] = VMT_CalHeapPos(RightNormE, RightComp, CalMethod);

BestE = [];
OKE = [];
LastBestE = zeros(StepSum);
IdealDisDiff = 0.025;

CompSum = abs(CompSide);
for i = 2: StepSum
    CompSum(i) = CompSum(i - 1) + CompSum(i);
end
U_0 = [0, (1: StepSum) - CompSum * OutputH * 2];    % 各个零势能点的位置

LoopTimes = -1;

if (~exist('BeginE', 'var') || size(BeginE, 1) == 0)
    ChangeRange = [ones(StepSum - 1, 1), 0.4 * ones(StepSum - 1, 1), RightNormE(2:end)' + 1];
else
    LoopTimes = 0;
    ChangeRange = [BeginE(2: StepSum)' - 0.3 * ones(StepSum - 1, 1), 0.1 * ones(StepSum - 1, 1), BeginE(2: StepSum)' + 0.3 * ones(StepSum - 1, 1)];
end

LoopAllRange = [];
NextBeginE = [];
MinMaxForceDiff = 5;

while(true)
    LoopTimes = LoopTimes + 1;
    fprintf('\n第%d次迭代求解：\n', LoopTimes);
    
    BestMinDisDiff = 0;
    OKSum = 0;
    CalSum = 0;

    AllCombination = GetSortedCombination(ChangeRange, LoopTimes ~= 0);
    ChangeSum = size(AllCombination, 1);
    fprintf('变化总数: %d\n', ChangeSum);

    
    if (LoopTimes ~= 0)
        if (mod(LoopTimes, 10) == 1)
            LoopAllRange = [LoopAllRange, zeros(StepSum - 1, 3 * 10)];
        end
        LoopAllRange(:, LoopTimes * 3 - 2: LoopTimes * 3) = ChangeRange;
    end


    t0 = tic;
    for ChangeNo = 1: ChangeSum
        ThisCombination = AllCombination(ChangeNo, :);
        LeftNormE = [1 ThisCombination];
        if (LoopTimes == 0)
            LeftNormE = LeftNormE + 0.2 * (0 : StepSum - 1);
        end
        AlreadyCalFlag = false;
        for i = 1: LoopTimes - 1
            if (all(ThisCombination' <= LoopAllRange(:, i * 3) + 10e-6) && all(ThisCombination' >= LoopAllRange(:, i * 3 - 2) - 10e-6))
                AlreadyCalFlag = true;
                break;
            end
        end
        if (AlreadyCalFlag)
            continue;
        end
        CalSum = CalSum + 1;
        [~, HeapPos_Left] = VMT_CalHeapPos(LeftNormE, LeftComp, CalMethod);
        LeftRightDiff = HeapPos_Left - HeapPos_Right + IfChanged * 4 * OutputH * (OriginStatus * 2 - 1);
        IfOK1 = ~any(((LeftRightDiff > 0) * 2 - 1) + GoalSequence * 2 - 1);
        % IfOK2 = true;
        IfOK2 = ~any((sign(LeftNormE - RightNormE) + GoalSequence) > 1);
        
       
        if (IfOK1 && IfOK2)
            
            ThisMaxForceDiff = 0;
            ChangeInfo = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)];
            Delta_HeapPos = ~[OriginStatus, GoalSequence(1: StepSum - 1)] * 2 * OutputH;
            Judge_HeapPos_L = HeapPos_Left + Delta_HeapPos;
            Judge_HeapPos_R = HeapPos_Right - Delta_HeapPos;
            for i = 1: StepSum
                if (ChangeInfo(i) == 0)
                    continue;
                end
                if (GoalSequence(i) == 1)
                    % 防止出现预先突跳
                    ThisDis = ((Judge_HeapPos_L(i) - U_0(i)) / (Judge_HeapPos_R(i) - U_0(i)) * RightNormE(i) - LeftNormE(i)) * Fm / OutputFm;
                else
                    ThisDis = ((Judge_HeapPos_R(i) - U_0(i)) / (Judge_HeapPos_L(i) - U_0(i)) * LeftNormE(i) - RightNormE(i)) * Fm / OutputFm;
                end
                ThisMaxForceDiff = max(ThisDis, ThisMaxForceDiff);
            end
    
            IfOK3 = ThisMaxForceDiff < 1;
            if (IfOK3)
                StepDisDiff = abs((HeapPos_Right - HeapPos_Left - 0.25 * (CompSide == 1)) .* CompSide);
                ThisMinDisDiff = min((StepDisDiff == 0) * 5 + StepDisDiff);
    
                OKSum = OKSum + 1;
                OKE = LeftNormE;
                if (ThisMaxForceDiff < MinMaxForceDiff && ThisMinDisDiff > IdealDisDiff) %  
                    BestE = LeftNormE;
                    MinMaxForceDiff = ThisMaxForceDiff;
                    BestMinDisDiff = ThisMinDisDiff;
                end
            end
        end
    end
    EpochTime = toc(t0);
    fprintf('实际计算数： %d\n', CalSum);
    fprintf('满足基本要求的设计数： %d\n', OKSum);
    fprintf('迭代总用时： %f  s\n', double(EpochTime));

    % 输出
    if (size(BestE))
        % 找到了最优的设计
        fprintf('最优设计：\n');
        for i = 1: StepSum
            fprintf('K%d = %.1f, ', i, BestE(i));
        end
        fprintf('\n');
        fprintf('最大力差异：%f\n', MinMaxForceDiff);
        NextBeginE = BestE;
        [BestRealE_Left, ~] = VMT_CalHeapPos(BestE, LeftComp, CalMethod);
    else
        % 没找到最优的设计
        fprintf('不存在较优的设计。\n');
        if (OKSum == 0)
            BestRealE_Left = zeros(1, StepSum);
            BestE = zeros(1, StepSum);
            break;
        else
            % 存在可行的设计
            fprintf('可行设计：\n');
            for i = 1: StepSum
                fprintf('K%d = %.1f, ', i, OKE(i));
            end
            fprintf('\n');
            NextBeginE = OKE;
            ChangeRange = [OKE(2: StepSum)' - 0.3 * ones(StepSum - 1, 1), 0.1 * ones(StepSum - 1, 1), OKE(2: StepSum)' + 0.3 * ones(StepSum - 1, 1)];
            continue;
        end
    end

    if (all(LastBestE == BestE))
        fprintf('此为最终最优设计。\n');
        break;
    end

    LastBestE = BestE;
    ChangeRange = [ NextBeginE(2: StepSum)' - 0.3 * ones(StepSum - 1, 1), ...
                    0.1 * ones(StepSum - 1, 1), ...
                    NextBeginE(2: StepSum)' + 0.3 * ones(StepSum - 1, 1)];
end
end