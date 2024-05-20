function GoalFunc = VMT_g_static(NormE, U_0, X_m, GoalSequence, OriginStatus, MaxNormE)
% Fmincon要优化的目标函数
    global g_CallTimes g_RunTime;
    persistent GoalFunc_static;
    g_CallTimes = g_CallTimes + 1;
    ThisRunTime = tic;
    StepSum = size(X_m, 2);
    TempNormE = sym('TempNormE_', [1, 2 * StepSum]);
    
    if (g_CallTimes == 1)
        % 代值，获得串联单元的峰值位置
        X_mL = X_m(1, :);
        X_mR = X_m(2, :);
        % X_mL = double(subs(X_m(1, :), TempNormE, NormE));
        % X_mR = double(subs(X_m(2, :), TempNormE, NormE));
    
        % Judge_X_m: 每一步用来判断哪侧先跳变，临时预计峰值位置。
        OutputH = 0.0625;
        Delta_HeapPos = ~[OriginStatus, GoalSequence(1: StepSum - 1)] * 2 * OutputH;
        Judge_X_mL = X_mL + Delta_HeapPos;
        Judge_X_mR = X_mR - Delta_HeapPos;
    
        % 优化函数
        GoalFunc_static = sym(0);
        % 这一步输出单元是否发生了跳变。
        ChangeInfo = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)];
        % 力的差异与位移的差异权重之比，用来调整优化策略。
        
    
        beta = 0.75;
        for i = 1: StepSum
            g_DisDiff = atan((Judge_X_mL(i) - Judge_X_mR(i)) * (2 * GoalSequence(i) - 1)) * 2 / pi;
            if (ChangeInfo(i) == 0)
                g_ForceDiff = 0;
            else
                if (GoalSequence(i) == 1)
                    g_ForceDiff = beta * (max((Judge_X_mL(i) - U_0(i)) / (Judge_X_mR(i) - U_0(i)) * TempNormE(StepSum + i) - TempNormE(i), 0));
                else
                    g_ForceDiff = beta * (max((Judge_X_mR(i) - U_0(i)) / (Judge_X_mL(i) - U_0(i)) * TempNormE(i) - TempNormE(StepSum + i), 0));
                end
                %g_ForceDiff = g_ForceDiff;
            end
            GoalFunc_static = GoalFunc_static + g_DisDiff + g_ForceDiff;
        end
        H_0 = 0.5;
        Comp_H_0 = H_0 - 2 * OutputH;
        Fm = 2*(1-(1+H_0^2)^(-1/3))^(3/2);
        Fm_Comp = 2*(1-(1+Comp_H_0^2)^(-1/3))^(3/2);
    
        CompSide = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)]; 
        LeftComp = CompSide == -1;
        RightComp = CompSide == 1;
        
        RealE_Left = TempNormE(1: StepSum) .* (1 + (LeftComp * (Fm/Fm_Comp - 1)));
        RealE_Right = TempNormE(StepSum + 1: 2 * StepSum) .* (1 + (RightComp * (Fm/Fm_Comp - 1)));
        g_FinalDisDiff = max((VMT_ConnectedGetU(RealE_Left, H_0 - LeftComp * 2 * OutputH, MaxNormE * Fm, ones(1, StepSum), 2)...
                        - VMT_ConnectedGetU(RealE_Right, H_0 - RightComp * 2 * OutputH, MaxNormE * Fm, ones(1, StepSum), 2)) * (1 - 2 * GoalSequence(StepSum)), 0);
        GoalFunc_static = GoalFunc_static + g_FinalDisDiff * StepSum;
    end

    GoalFunc = double(subs(GoalFunc_static, TempNormE, NormE));
    g_RunTime = g_RunTime + double(toc(ThisRunTime));
end