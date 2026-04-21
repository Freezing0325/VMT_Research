function VMT_ReportConfig(fileID, BestE, GoalSequence, OriginStatus, CalMethod)
    StepSum = size(GoalSequence, 2);
    if (size(GoalSequence,1) == 1)
        GoalSequence_Hat = [OriginStatus, GoalSequence(1: StepSum - 1)];
        % 需要施加补偿的一侧，0：不需要，-1：左侧，1：右侧。
        CompSide = GoalSequence - GoalSequence_Hat; 
        LeftComp = CompSide == -1;
        RightComp = CompSide == 1;
    else
        GoalSequence_Hat = [OriginStatus, GoalSequence(1, 1: StepSum - 1)];
        LeftComp = GoalSequence(2, :);
        RightComp = GoalSequence(3, :);
        GoalSequence = GoalSequence(1, :);
        CompSide = GoalSequence - GoalSequence_Hat; 
    end
    global Output_h a Normal_h;
    if (isempty(Output_h))
        VMT_Init();
    end
    H_0 = Normal_h / a;
    OutputH = Output_h / a;
    
    if (size(BestE, 2) == 2 * (StepSum-1) && size(BestE, 1) == 1)
        LeftNormE = [1,BestE(1: StepSum - 1)];
        RightNormE = [1,BestE(StepSum: 2 * (StepSum - 1))];
    else
        TempBestE = BestE';
        LeftNormE = TempBestE(1: StepSum)';
        RightNormE = TempBestE(StepSum+1: 2*StepSum)';
    end
    BestNormE = [LeftNormE, RightNormE];
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
    


    fprintf(fileID, '归一化刚度：\n');
    for i = 1: 2 * StepSum
        fprintf(fileID, '%.4f  ', BestNormE(i));
    end
    fprintf(fileID, '\n');
    fprintf(fileID, '实际刚度：\n');
    for i = 1: 2 * StepSum
        fprintf(fileID, '%.4f  ', RealE(floor((i-1)/StepSum)+1, mod(i-1,StepSum)+1));
    end
    fprintf(fileID, '\n');

    [PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE , LeftComp, RightComp, OriginStatus, CalMethod, []);
    fprintf(fileID, '预期序列：\n');
    for i = 1: StepSum
        fprintf(fileID, '%d  ', PredSequence(i));
    end
    fprintf(fileID, '\n');
    fprintf(fileID, '最大力差异：%f\n', MaxForceDiff);

    [MinDispDiff, MDDIndex] = min(abs(Judge_H(1,2:end)-Judge_H(2,2:end)));
    fprintf(fileID, '最小位移差异：%f，出现在第%d阶段\n', MinDispDiff, MDDIndex + 1);
    
    [Fm, ~] = VMT_SingleGetFm(1, H_0, CalMethod);
    FinalDisDiff = (VMT_ConnectedGetU(RealE(1,:), H_0 - LeftComp * 2 * OutputH, BestMaxNormE * Fm, ones(1, StepSum), CalMethod)...
                            - VMT_ConnectedGetU(RealE(2,:), H_0 - RightComp * 2 * OutputH, BestMaxNormE * Fm, ones(1, StepSum), CalMethod)) * (1 - 2 * GoalSequence(StepSum));
    fprintf(fileID, '最终位移差异：%f\n', FinalDisDiff/OutputH); % 这个位移差异是考虑到最终状态时的结果，负值更稳定。
    
end