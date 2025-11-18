function [g, h] = VMT_con_static(NormE, X_m, GoalSequence, OriginStatus, U_0, CalMethod, MinDisDiff, MinStepWall, MaxNormE, MaxOutDisDiff, MaxFDiff)
% Fmincon要满足的约束条件，g<=0,h=0
    global con_CallTimes con_RunTime;
    persistent g_static h_static;
    
    ThisRunTime = tic;
    StepSum = size(X_m, 2);
    TempNormE = sym('TempNormE_', [1, 2 * StepSum]);
    
    if (con_CallTimes == 0)

        % 代值，获得串联单元的峰值位置
        
        X_mL = X_m(1, :);
        X_mR = X_m(2, :);

        global a Normal_h OutputEA_ka Output_h Output_a;
        if (isempty(Output_h))
            VMT_Init();
        end
        OutputH = Output_h / a;
        OutputA = Output_a / a;
        [OutputFm, OutputHm] = VMT_SingleGetFm(OutputEA_ka, OutputH/OutputA, CalMethod);
        H_0 = Normal_h / a;
        Comp_H_0 = H_0 - 2 * OutputH;
        
        [Fm, ~] = VMT_SingleGetFm(1, H_0, CalMethod);
        [Fm_Comp, ~] = VMT_SingleGetFm(1, Comp_H_0, CalMethod);
        GoalSequence_Hat = [OriginStatus, GoalSequence(1: StepSum - 1)];
        % 这一步输出单元是否发生了跳变。1:0→1,-1:1→0
        ChangeInfo = GoalSequence - GoalSequence_Hat;
    
        % Judge_X_m: 每一步用来判断哪侧先跳变，临时预计峰值位置。实际的位置可能与此不同，因为输出单元的效应。
        Delta_HeapPos = (OriginStatus - GoalSequence_Hat) * 2 * OutputH;
        Judge_X_mL = X_mL + Delta_HeapPos;
        Judge_X_mR = X_mR - Delta_HeapPos;
    
        % Real_X_m: 在序列已经确定的前提下，两侧真实的峰值位置，也就是曲线图上的实际位置。
        Real_X_mL = Judge_X_mL + (ChangeInfo == -1) * 2 * OutputH;
        Real_X_mR = Judge_X_mR + (ChangeInfo == 1) * 2 * OutputH;
        
    
        % 第一个约束，对位移的约束：要求峰值位置和想要设计的变形序列相匹配，同时留下一定的容许误差空间MinDisDiff。
        g_DisDiff = ((Judge_X_mL - Judge_X_mR) .* (GoalSequence * 2 - 1) + [0, MinDisDiff * ones(1, StepSum - 1)]);
    
        % 第二个约束，对相邻两个bit之间的位置的约束：要求每一个bit严格只有一个左侧和右侧单元发生突跳。
        % 序列发生变化的时候，不会有这方面的影响，但序列不变的时候，例如说要"保持"为1，其它的性质满足的条件是：
        % 1.上一个bit，X_mL(i - 1) < X_mR(i - 1)
        % 2.这一个bit，X_mL(i) < X_mR(i)
        % 3.隐含性质，X_mL(i - 1) < X_mL(i); X_mR(i - 1) < X_mR(i)
        % 这个过程没有对  X_mR(i - 1)  和  X_mL(i)  做出约束， 
        % 因而如果 X_mL(i - 1) < X_mL(i) < X_mR(i - 1) < X_mR(i)，那么不同的bit之间就会发生相互干扰。
        % 所以强制要求，Judge_X_mL(i) - Real_X_mR(i - 1) > MinStepWall，从而使得不同bit之前产生"壁垒"。
    
        g_StepWall = sym(zeros(1, StepSum - 1));
    
        % 第三个约束，对突跳前力差异的约束：如果这一步输出序列要发生一次转变，那么转变必须严格发生在一侧单元先行突跳的位置。
        % 如果在两侧单元都处于"爬坡"阶段的时候就已经产生了较大的力的差异，那么就有可能提前发生转变，不符合设计要求。
     
        AllForceDiff = sym(zeros(1, StepSum));
    
        
        for i = 1: StepSum
            if (i > 1 && ChangeInfo(i) == 0)
                if (GoalSequence(i) == 1)
                    g_StepWall(i - 1) = Real_X_mR(i - 1) - Judge_X_mL(i) + MinStepWall;
                else
                    g_StepWall(i - 1) = Real_X_mL(i - 1) - Judge_X_mR(i) + MinStepWall;
                end
            end
            k1 = TempNormE(i) / (Judge_X_mL(i) - U_0(i)) * Fm;
            k2 = TempNormE(StepSum + i) / (Judge_X_mR(i) - U_0(i)) * Fm;
            HmOutputFm = (k1 + k2) * (OutputH - OutputHm) + OutputFm;
            if (GoalSequence(i) == 1)
                ThisDis = ((Judge_X_mL(i) - U_0(i)) * k2 - TempNormE(i) * Fm) / HmOutputFm;
            else
                ThisDis = (TempNormE(StepSum + i) * Fm - (Judge_X_mR(i) - U_0(i)) * k1) / HmOutputFm;
            end
            ThisDis = ThisDis * (1 - 2 * GoalSequence_Hat(i));
            AllForceDiff(i) = ThisDis - MaxFDiff;
        end
    
        g_ForceDiff = AllForceDiff;
        g_StepWall = g_StepWall(g_StepWall ~= 0);
    
        
    
        CompSide = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)]; 
        LeftComp = CompSide == -1;
        RightComp = CompSide == 1;
        RealE_Left = TempNormE(1: StepSum) .* (1 + (LeftComp * (Fm/Fm_Comp - 1)));
        RealE_Right = TempNormE(StepSum + 1: 2 * StepSum) .* (1 + (RightComp * (Fm/Fm_Comp - 1)));
        g_FinalDisDiff = (VMT_ConnectedGetU(RealE_Left, H_0 - LeftComp * 2 * OutputH, MaxNormE * Fm, ones(1, StepSum), 2)...
                        - VMT_ConnectedGetU(RealE_Right, H_0 - RightComp * 2 * OutputH, MaxNormE * Fm, ones(1, StepSum), 2)) * (1 - 2 * GoalSequence(StepSum)) - MaxOutDisDiff;
    
        
    
        % MinNearETimes = 0.98;
        % g_NearETimes = [NormE(1: StepSum - 1) ./ NormE(2: StepSum), NormE(StepSum + 1: 2 * StepSum - 1) ./ NormE(StepSum + 2: 2 * StepSum)] - MinNearETimes;
        % g_StepWall.'; 
        g_static = [g_DisDiff.'; g_ForceDiff.'; g_FinalDisDiff];  %; g_NearETimes.'
        h_static = [];
    end

    g = double(subs(g_static, TempNormE, NormE));
    h = double(subs(h_static, TempNormE, NormE));
    con_CallTimes = con_CallTimes + 1;
    con_RunTime = con_RunTime + double(toc(ThisRunTime));
    
end