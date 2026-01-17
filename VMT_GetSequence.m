function [PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, ActiveStatus)
%VMT_GetSequence       预测当前配置的VMT所对应的变化序列。
%
%   输入参数：
%   LeftNormE       左侧所有单元的归一化刚度，为EA=kl0的值
%   RightNormE      右侧所有单元的归一化刚度，为EA=kl0的值
%   LeftComp        左侧补偿情况
%   RightComp       右侧补偿情况
%   OriginStatus    初始状态
%   CalMethod       计算算法，1：求解非线性方程；2：线性模型：泰勒级数2阶拟合；3：非线性模型：泰勒级数3阶拟合。
%   ActiveStatus
%   两侧单元激活情况，第一行为左侧，第二行为右侧。不输入时，认为所有单元均激活。也可以只输入不激活的单元序号，分为两行。
%
%   输出：[PredSequence, MaxForceDiff]
%   PredSequence    预测的序列
%   如果出现-2ab，则说明在该位置状态将变成a，但会出现突跳顺序或是补偿机制不正确的现象，尽管突跳顺序决定为状态a，但几何形状使得状态b更稳定。
%   如果出现-3ab，则说明在该位置由于力的差异大于输出单元的突跳阈值，会发生意外的突跳，突跳为a，但随后的突跳顺序会产生b状态。
%   MaxForceDiff    会影响序列变化的最大的力差异，用输出单元的突跳最大力得到比值

    SingleSideComp = true;  % 仅进行单侧的补偿
    global a Normal_h OutputEA_ka Output_h Output_a;
    if (isempty(Output_h))
        VMT_Init();
    end
    OutputH = Output_h / a;
    OutputA = Output_a / a;
    [OutputFm, OutputHm] = VMT_SingleGetFm(OutputEA_ka, OutputH / OutputA, CalMethod);
    [Fm, ~] = VMT_SingleGetFm(1, Normal_h / a, CalMethod);

    UnitSum = size(LeftNormE, 2);
    if (~exist('ActiveStatus', 'var') || size(ActiveStatus, 1) == 0)
        ActiveStatus = ones(2, UnitSum);
    end
    if (size(ActiveStatus, 2) < UnitSum)
        TempActiveStatus = ones(2, UnitSum);
        TempActiveStatus(1, ActiveStatus(1,:)) = 0;
        TempActiveStatus(2, ActiveStatus(2,:)) = 0;
        ActiveStatus = TempActiveStatus;
    end

    LeftActiveStatus = ActiveStatus(1, :);
    RightActiveStatus = ActiveStatus(2, :);
    LeftActiveIndex = find(LeftActiveStatus ~= 0);
    RightActiveIndex = find(RightActiveStatus ~= 0);
    LeftStepSum = size(LeftActiveIndex, 2);
    RightStepSum = size(RightActiveIndex, 2);
    ActiveCompInfo = [LeftComp(LeftActiveIndex);RightComp(RightActiveIndex)];
    if (LeftStepSum ~= RightStepSum)
        error('左右两侧激活数不同！\n');
    end
    StepSum = LeftStepSum;
    
    % 串联结构本身压缩时的势能零点位置
    MaterialU_0 = zeros(2, StepSum+1);
    [~, HeapPos_Left] = VMT_CalHeapPos(LeftNormE, LeftComp, CalMethod, LeftActiveStatus);
    [~, HeapPos_Right] = VMT_CalHeapPos(RightNormE, RightComp, CalMethod, RightActiveStatus);
    NowStatus = OriginStatus;
    GeometryBalance = OriginStatus;
    GeometryStatus = OriginStatus;
    PredSequence = zeros(1, StepSum);
    % 在某阶段前的变形序列发生之后，这一步用来判断的峰值位置
    U_0 = zeros(2, StepSum+1);
    Judge_HeapPos_L = HeapPos_Left;
    Judge_HeapPos_R = HeapPos_Right;
    MaxForceDiff = 0;
    for i = 1: StepSum
        MaterialU_0(:,i+1) = MaterialU_0(:,i) + 2 * (Normal_h / a - ActiveCompInfo(:,i) * OutputH * (1 + SingleSideComp));
        % 突跳发生前，先看是否稳定，此时达到的最大压力差是否支持维持在这一状态
        ThisDeltaHeapPos = (OriginStatus - NowStatus) * 2 * OutputH;
        Judge_HeapPos_L(i) = Judge_HeapPos_L(i) + ThisDeltaHeapPos;
        Judge_HeapPos_R(i) = Judge_HeapPos_R(i) - ThisDeltaHeapPos;
        PredSequence(i) = Judge_HeapPos_L(i) < Judge_HeapPos_R(i);
        % PredSequence(i)表征了突跳发生的顺序。
        % 计算的压力差值是右侧力-左侧力
        k1 = LeftNormE(LeftActiveIndex(i)) / (Judge_HeapPos_L(i) - U_0(1, i)) * Fm;
        k2 = RightNormE(RightActiveIndex(i)) / (Judge_HeapPos_R(i) - U_0(2, i)) * Fm;
        % ZeroOutputFm = (k1 + k2) * OutputH;
        % 这种计算方式是要求变形到输出单元完全呈一条线（不稳定平衡点）时，力的大小仍然要大于零，才能推动其切换
        HmOutputFm = (k1 + k2) * (OutputH - OutputHm) + OutputFm;
        % 这种计算方式是要求变形到输出单元达到其切换的极大值点时，力的大小仍然要大于阈值力，才能推动其切换
        % 在 HmOutputFm > ZeroOutputFm 时，HmOutputFm 最为可靠，是最合理的判据。
        % 但ZeroOutputFm > HmOutputFm时就不好说了。
        % fprintf('三种计算方式的结果：\n%f\n%f\n%f\n', OutputFm, ZeroOutputFm, HmOutputFm);

        
        if (PredSequence(i) == 1)   % 本次突跳左侧早于右侧
            %ThisDis = (Judge_HeapPos_L(i) - U_0(2, i) - (OutputH - OutputHm)) * (k2 - k1) / HmOutputFm;
            ThisDis = ((Judge_HeapPos_L(i) - U_0(2, i)) * k2 - LeftNormE(LeftActiveIndex(i)) * Fm) / HmOutputFm;
        else % 本次突跳右侧早于左侧
            %ThisDis = (Judge_HeapPos_R(i) - U_0(2, i) - (OutputH - OutputHm)) * (k2 - k1) / HmOutputFm;
            ThisDis = (RightNormE(RightActiveIndex(i)) * Fm - (Judge_HeapPos_R(i) - U_0(1, i)) * k1) / HmOutputFm;
        end
        % 当前状态为0时，压力差为负值安全；当前状态为1时，压力差为正值安全
        ThisDis = ThisDis * (1 - 2 * NowStatus);
        % 处理之后，ThisDis为负值安全，为正值不安全
        MaxForceDiff = max(ThisDis, MaxForceDiff);
        % 不安全，将提前切换至1-NowStatus状态。
        if (ThisDis >= 1)
            NowStatus = 1-NowStatus;
            PredSequence(i) = -300 - NowStatus*10;
            ThisDeltaHeapPos = (OriginStatus - NowStatus) * 2 * OutputH;
            Judge_HeapPos_L(i) = Judge_HeapPos_L(i) + ThisDeltaHeapPos;
            Judge_HeapPos_R(i) = Judge_HeapPos_R(i) - ThisDeltaHeapPos;
            PredSequence(i) = PredSequence(i) - (Judge_HeapPos_L(i) < Judge_HeapPos_R(i));
        end
        
        GeometryBalance = GeometryBalance + RightComp(RightActiveIndex(i)) - LeftComp(LeftActiveIndex(i));
        GeometryStatus = min(max(GeometryBalance, 0), 1);
        % 如果根据突跳顺序预测的状态与根据补偿单元预测的状态不一致
        ThisPredSeq = mod(abs(PredSequence(i)),10);
        if (ThisPredSeq ~= GeometryStatus)
            if (PredSequence(i) <= -300)
                PredSequence(i) = PredSequence(i) * 1000 - (200 + ThisPredSeq * 10 + GeometryStatus);
            elseif (i > 1)
                PredSequence(i) = -(200 + ThisPredSeq * 10 + GeometryStatus);
            else
                PredSequence(i) = GeometryStatus;
                ThisPredSeq = GeometryStatus;
            end
        end
        NowStatus = ThisPredSeq;
        U_0(:,i+1) = MaterialU_0(:,i+1) + [-1; 1] * (NowStatus-OriginStatus) * OutputH * 2;
    end
    UnstableSwitchIndex = find(PredSequence < 0);
    if (size(UnstableSwitchIndex,2) > 0 && UnstableSwitchIndex(1) ~= StepSum)
        MaxForceDiff = max(MaxForceDiff, 1);
        return
    end
end
