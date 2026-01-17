function GoalFunc = VMT_g_static(FullNormE, GoalSequence, OriginStatus, CalMethod, Optimizer)
% Fmincon要优化的目标函数
    global g_CallTimes g_RunTime;
    persistent GoalFunc_static;
    
    ThisRunTime = tic;
    StepSum = size(GoalSequence, 2);
    TempNormE = sym('TempNormE_', [1, 2 * StepSum]);

    if (g_CallTimes == 0)
        GoalSequence_Hat = [OriginStatus, GoalSequence(1: StepSum - 1)];
        CompSide = GoalSequence - GoalSequence_Hat; 
        LeftComp = CompSide == -1;
        RightComp = CompSide == 1;
        CompSum = abs(CompSide);
        for i = 2: StepSum
            CompSum(i) = CompSum(i - 1) + CompSum(i);
        end
        % 用符号表示的各个峰值位置，是归一化刚度的函数
        [~, X_mL] = VMT_CalHeapPos_2(TempNormE(1: StepSum), LeftComp);
        [~, X_mR] = VMT_CalHeapPos_2(TempNormE(StepSum + 1: 2 * StepSum), RightComp);
        % X_m = [X_mL; X_mR];
        % X_mL = X_m(1, :);
        % X_mR = X_m(2, :);
        % X_mL = double(subs(X_m(1, :), TempNormE, NormE));
        % X_mR = double(subs(X_m(2, :), TempNormE, NormE));
        global a Normal_h Output_h OutputEA_ka Output_a;
        if (isempty(Output_h))
            VMT_Init();
        end
        OutputH = Output_h / a;
        OutputA = Output_a / a;
        [OutputFm, ~] = VMT_SingleGetFm(OutputEA_ka, OutputH/OutputA, CalMethod);
        H_0 = Normal_h / a;
        U_0 = [0, 2*(1: StepSum)*H_0 - CompSum * OutputH * 2];
        Comp_H_0 = H_0 - 2 * OutputH;
        [Fm, ~] = VMT_SingleGetFm(1, H_0, CalMethod);
        [Fm_Comp, ~] = VMT_SingleGetFm(1, Comp_H_0, CalMethod);
        GoalSequence_Hat = [OriginStatus, GoalSequence(1: StepSum - 1)];
        % 这一步输出单元是否发生了跳变。
        ChangeInfo = GoalSequence - GoalSequence_Hat;
    
        % Judge_X_m: 每一步用来判断哪侧先跳变，临时预计峰值位置。
        Delta_HeapPos = (OriginStatus - GoalSequence_Hat) * 2 * OutputH;
        Judge_X_mL = X_mL + Delta_HeapPos;
        Judge_X_mR = X_mR - Delta_HeapPos;
    
        % 优化函数
        GoalFunc_static = sym(0);
        
        % 力的差异与位移的差异权重之比，用来调整优化策略。
        
    
        beta = 0.2 / OutputFm;
        for i = 1: StepSum
            g_DisDiff = atan((Judge_X_mL(i) - Judge_X_mR(i)) * (2 * GoalSequence(i) - 1)) * 2 / pi;
            if (GoalSequence(i) == 1)
                ThisDis = ((Judge_X_mL(i) - U_0(i)) * TempNormE(StepSum + i) / (Judge_X_mR(i) - U_0(i)) - TempNormE(i));
            else
                ThisDis = (TempNormE(StepSum + i) - (Judge_X_mR(i) - U_0(i)) * TempNormE(i) / (Judge_X_mL(i) - U_0(i)));
            end
            g_ForceDiff = beta * max(ThisDis * (1 - 2 * GoalSequence_Hat(i)), 0);
            GoalFunc_static = GoalFunc_static + g_DisDiff + g_ForceDiff;
        end
        
        CompSide = ChangeInfo; 
        LeftComp = CompSide == -1;
        RightComp = CompSide == 1;
        
        RealE_Left = TempNormE(1: StepSum) .* (1 + (LeftComp * (Fm/Fm_Comp - 1)));
        RealE_Right = TempNormE(StepSum + 1: 2 * StepSum) .* (1 + (RightComp * (Fm/Fm_Comp - 1)));
        g_FinalDisDiff = max((VMT_ConnectedGetU(RealE_Left, H_0 - LeftComp * 2 * OutputH, Optimizer.MaxNormE * Fm, ones(1, StepSum), 2)...
                        - VMT_ConnectedGetU(RealE_Right, H_0 - RightComp * 2 * OutputH, Optimizer.MaxNormE * Fm, ones(1, StepSum), 2)) * (1 - 2 * GoalSequence(StepSum)), 0);
        GoalFunc_static = GoalFunc_static + g_FinalDisDiff * StepSum;
    end
    
    GoalFunc = double(subs(GoalFunc_static, TempNormE, FullNormE));
    g_CallTimes = g_CallTimes + 1;
    g_RunTime = g_RunTime + double(toc(ThisRunTime));
end