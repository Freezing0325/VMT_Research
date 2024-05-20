RightSideNormK = [1 1.5 3 3.5];
OriginStatus = 1;   
GoalSequence = [1 0 1 0];   % 1：左高右低，0：左低右高。
StepSum = size(GoalSequence, 2);    % 序列长度
OutputH = 0.0625;

IfChanged = [0, GoalSequence(1: end - 1) ~= OriginStatus];  % 目前的状态是否不同于初始状态
CompensationSide = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)]; % 需要施加补偿的一侧，0：不需要，-1：左侧，1：右侧。
LeftComp = CompensationSide == -1;
RightComp = CompensationSide == 1;
CompSum = abs(CompensationSide);
for i = 2: StepSum
    CompSum(i) = CompSum(i - 1) + CompSum(i);
end
U_0 = [0, (1: StepSum) - CompSum * OutputH * 2];    % 各个零势能点的位置


[RealK_Right, HeapPos_Right] = CalHeapPos(RightSideNormK, RightComp);
% LeftSideNormK = [1 2 2.5 4];
% [ActualRealK_Left, ActualHeapPos_Left] = CalHeapPos(LeftSideNormK, LeftComp);


%%
BestK = [];
Rob = 0.025;
BestRob = 0;
MinDis = 5;
OKSum = 0;
for K2 = 1.1: 0.1: 4
    for K3 = K2 + 0.1: 0.1: 6
        for K4 = K3 + 0.1: 0.1: 8
            LeftSideNormK = [1 K2 K3 K4];
            [RealK_Left, HeapPos_Left] = CalHeapPos(LeftSideNormK, LeftComp);
            LeftRightDiff = HeapPos_Left - HeapPos_Right + IfChanged * 4 * OutputH * (OriginStatus * 2 - 1);
            IfOK1 = ~any(((LeftRightDiff > 0) * 2 - 1) + GoalSequence * 2 - 1);
            IfOK2 = ~any((sign(LeftSideNormK - RightSideNormK) + GoalSequence) > 1);
            % if (all(LeftSideNormK == RightSideNormK))
            %     fprintf('yeah\n');
            % end
            if (IfOK1 && IfOK2)
                StepRob = abs((HeapPos_Right - HeapPos_Left - 0.25 * (CompensationSide == 1)) .* CompensationSide);
                ThisRob = min((StepRob == 0) * 5 + StepRob);

                ThisDis = ((LeftSideNormK .* (HeapPos_Right - U_0(1: end - 1))) ./ (RightSideNormK .* (HeapPos_Left - U_0(1: end - 1)))) .^ (1 - 2 * GoalSequence(end)) - 1;
                AllNormK = [LeftSideNormK; RightSideNormK];
                ThisDis = ThisDis .* AllNormK([2 - GoalSequence + (0 : StepSum - 1) * 2]);
                ThisMaxDis = max(ThisDis);
                OKSum = OKSum + 1;

                if (ThisMaxDis < MinDis && ThisRob > Rob)
                    BestK = LeftSideNormK;
                    MinDis = ThisMaxDis;
                    BestRob = ThisRob;
                end
                % fprintf('K1 = %d, K2 = %d, K3 = %d, K4 = %d, OK!\n', LeftSideNormK(1), LeftSideNormK(2), LeftSideNormK(3), LeftSideNormK(4));
            end
        end
    end
end
fprintf('OKSum = %d\n', OKSum);
if (size(BestK))
    fprintf('K1 = %d, K2 = %d, K3 = %d, K4 = %d, Best!\n', BestK(1), BestK(2), BestK(3), BestK(4));
    [BestRealK_Left, BestHeapPos_Left] = CalHeapPos(BestK, LeftComp);
else
    fprintf('No Best Solution.\n');
end


%%
function RealHeapPos = GetRealHeapPos(SingleHeapPos, SideChange)
    % SideChange：这一侧受到拉伸为1，压缩为-1，不变为0.
    OutputH = 0.0625;
    RealHeapPos = SingleHeapPos;
    StepSum = size(SideChange, 2);
    for i = 1: StepSum
        if (SideChange(i) ~= 0)
            RealHeapPos(i + (SideChange(i) == -1): end) = RealHeapPos(i + (SideChange(i) == -1): end) + SideChange(i) * 2 * OutputH;
        end
    end
end

function [RealK, HeapPos] = CalHeapPos(NormK, CompensationStatus)
    CompensationH_0 = 0.375;
    H_0 = 0.5;
    HMat = [H_0, CompensationH_0];
    % theta_0_Mat = atan(HMat);
    % theta_m_Mat = atan(HMat ./ sqrt(2*HMat.^2 + 3));
    U_snap_Mat = HMat - HMat ./ sqrt(2*HMat.^2 + 3);
    Func_Mat = zeros(2, 3);
    Func_Mat(:, 1) = 2*HMat.^2 ./ (HMat.^2 + 1).^(3/2);
    Func_Mat(:, 2) = (3*HMat.*(HMat.^2 - 1)) ./ (HMat.^2 + 1).^(5/2);
    Func_Mat(:, 3) = (4*HMat.^4 - 10*HMat.^2 + 1) ./ (HMat.^2 + 1).^(7/2);
    Func_snap_Mat = 2 * HMat.^3 ./ (3 * sqrt(3) .* sqrt(HMat.^2 + 1));
    % L_0_Mat = sqrt(HMat.^2 + 1);
    F_m_Mat = Func_snap_Mat;
    
    CorrK = F_m_Mat(1) / F_m_Mat(2);
    RealK = (CorrK - 1) * CompensationStatus .* NormK + NormK;

    [~, SortIndex] = sort(NormK);
    UnitSum = size(NormK, 2);
    HeapPos = zeros(1, UnitSum);
    for i = 1: UnitSum
        IfThisComp = CompensationStatus(SortIndex(i));
        U_all = zeros(1, UnitSum);
        for j = 1: i - 1
            IfThatComp = CompensationStatus(SortIndex(j));
            Eqn_a = -Func_Mat(IfThatComp + 1, 3);
            Eqn_b = -Func_Mat(IfThatComp + 1, 2);
            Eqn_c = -Func_Mat(IfThatComp + 1, 1);
            Eqn_d = -RealK(SortIndex(i)) / RealK(SortIndex(j)) * Func_snap_Mat(IfThisComp + 1);
            Eqn_p = Eqn_c / Eqn_a - (Eqn_b / Eqn_a)^2 / 3;
            Eqn_q = Eqn_d / Eqn_a + 2 * (Eqn_b / (3 * Eqn_a))^3  - Eqn_b * Eqn_c / (3 * Eqn_a^2);
            Eqn_r = sqrt(-(Eqn_p/3)^3);
            Eqn_theta = acos(-Eqn_q / (2 * Eqn_r)) / 3;
            U_all(j) = 2 * HMat(IfThatComp + 1) - 2 * Eqn_r^(1/3) * cos(Eqn_theta - 2 * pi / 3) + Eqn_b / (3 * Eqn_a);
        end
        U_all(i) = U_snap_Mat(IfThisComp + 1);
        for j = i + 1: UnitSum
            IfThatComp = CompensationStatus(SortIndex(j));
            Eqn_a = -Func_Mat(IfThatComp + 1, 3);
            Eqn_b = -Func_Mat(IfThatComp + 1, 2);
            Eqn_c = -Func_Mat(IfThatComp + 1, 1);
            Eqn_d = RealK(SortIndex(i)) / RealK(SortIndex(j)) * Func_snap_Mat(IfThisComp + 1);
            Eqn_p = Eqn_c / Eqn_a - (Eqn_b / Eqn_a)^2 / 3;
            Eqn_q = Eqn_d / Eqn_a + 2 * (Eqn_b / (3 * Eqn_a))^3  - Eqn_b * Eqn_c / (3 * Eqn_a^2);
            Eqn_r = sqrt(-(Eqn_p/3)^3);
            Eqn_theta = acos(-Eqn_q / (2 * Eqn_r)) / 3;
            U_all(j) = 2 * Eqn_r^(1/3) * cos(Eqn_theta - 2 * pi / 3) - Eqn_b / (3 * Eqn_a);
        end
        HeapPos(i) = sum(U_all);
    end
end



