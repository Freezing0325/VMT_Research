function [RealE, HeapPos] = VMT_CalHeapPos_2(NormE, CompStatus)
%VMT_CalHeapPos_2     用二、三阶泰勒级数、线性模型计算当前一侧配置的串联VMT所对应的峰值位置
%
%   NormK           所有单元的归一化刚度
%   CompStatus      补偿情况，1表示有补偿，0表示无补偿

    global a Normal_h CompPlus_h CompMinus_h;
    if (isempty(Normal_h))
        VMT_Init();
    end
    NormalH_0 = Normal_h / a;     %0.5
    CompMinusH_0 = CompMinus_h / a;  %0.375
    CompPlusH_0 = CompPlus_h / a;   %0.625
    
    HMat = [CompPlusH_0, NormalH_0, CompMinusH_0];
    
    [F_snap_Mat, U_snap_Mat] = VMT_SingleGetFm(1, HMat, 2);
    

    CorrE = F_snap_Mat(2) ./ F_snap_Mat;
    
    [NormE, Index] = sort(NormE);
    CompStatus = CompStatus(Index);

    RealE = NormE .* CorrE(CompStatus + 2);

    UnitSum = size(NormE, 2);
    HeapPos = sym(zeros(1, UnitSum));

    for i = 1: UnitSum
        IfThisComp = CompStatus(i);
        U_all = sym(zeros(1, UnitSum));
        for j = 1: UnitSum
            if (j == i)
                U_all(j) = U_snap_Mat(IfThisComp + 2);
                continue;
            end
            IfThatComp = CompStatus(j);
            IsDown = (j < i);
            % U_all_j_2 = VMT_SingleGetU(RealE(j), HMat(IfThatComp + 2), RealE(i) * F_snap_Mat(IfThisComp + 2), IsDown, 2);
            % if (IsDown)
            % U_all(j) = U_all_j_2;
            % else
            %     U_all_j_1 = VMT_SingleGetU(RealE(j), HMat(IfThatComp + 2), RealE(i) * F_snap_Mat(IfThisComp + 2), 0, -2);
            %     kEq_ratio = NormEA_ka(j) / NormEA_ka(i);
            %     alpha_1 = 0.5 * tanh(40*(kEq_ratio - 1.25)) + 0.5;
            %     U_all(j) = U_all_j_1 * (1-alpha_1) + U_all_j_2 * alpha_1;
            % end
            if (j == i + 1)
                U_all(j) = VMT_SingleGetU(RealE(j), HMat(IfThatComp + 2), RealE(i) * F_snap_Mat(IfThisComp + 2), 0, -2);
            else
                U_all(j) = VMT_SingleGetU(RealE(j), HMat(IfThatComp + 2), RealE(i) * F_snap_Mat(IfThisComp + 2), IsDown, 2);
            end
        end
        HeapPos(i) = sum(U_all);
    end
end