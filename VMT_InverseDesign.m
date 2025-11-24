function [BestE, Bestg] = VMT_InverseDesign(GoalSequence, OriginStatus, CalMethod, BeginNormE, Optimizer)

    % 目标序列
    % GoalSequence = [1 0 0 1 0 1 1 1 0];
    StepSum = size(GoalSequence, 2);
    
    
    
    % CalMethod = 2;
    
    % 迭代开始的归一化刚度配置
    % BeginNormE = [1.5 2, 1.5 2];
    % BeginNormE = [1.5     2       2.5     3         3.5     4       4.5     5, ...
    %               1.5     2       2.5     3         3.5     4       4.5     5];
    
    

    % 获得各个零势能点的位置
    global Output_h;
    if (isempty(Output_h))
        VMT_Init();
    end
   
    
    
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
    % [-M; -M; -M; ...; -M; -M; -M; -M; ...; -M];
    
    A_SortE_temp = [zeros(1, StepSum - 1); eye(StepSum - 2), zeros(StepSum - 2, 1)] - eye(StepSum - 1);
    A_SortE = [A_SortE_temp, zeros(StepSum - 1); zeros(StepSum - 1), A_SortE_temp];
    B_SortE = [-1; zeros(StepSum - 2, 1); -1; zeros(StepSum - 2, 1)] - Optimizer.MinNormEDiff * ones(2 * (StepSum - 1), 1);
    global g_CallTimes con_CallTimes g_RunTime con_RunTime;
    g_RunTime = 0;
    con_RunTime = 0;
    g_CallTimes = 0;
    con_CallTimes = 0;
    fprintf('开始迭代计算，过程可能要很久。\n')
    AllRunTime = tic;
    [BestE, Bestg] = fmincon(@(NormE)VMT_g_static([1, NormE(1: StepSum - 1), 1, NormE(StepSum: 2 * (StepSum - 1))], GoalSequence, OriginStatus, CalMethod, Optimizer), ...
                            BeginNormE, A_SortE, B_SortE, [], [], 1 * ones(1, 2*StepSum-2), Optimizer.MaxNormE * ones(1, 2*StepSum-2), ...
                            @(NormE)VMT_con_static([1, NormE(1: StepSum - 1), 1, NormE(StepSum: 2 * (StepSum - 1))], GoalSequence, OriginStatus, CalMethod, Optimizer));
    toc(AllRunTime);
    
    % QQ_Report('1603441246', 'Matlab算完了噢~');
end


