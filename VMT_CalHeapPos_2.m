function [RealE, HeapPos] = VMT_CalHeapPos_2(NormE, CompStatus)
%VMT_CalHeapPos_2     用二、三阶泰勒级数、线性模型计算当前一侧配置的串联VMT所对应的峰值位置
%
%   NormK           所有单元的归一化刚度
%   CompStatus      补偿情况，1表示有补偿，0表示无补偿

    CompensationH_0 = 0.375;
    NormalH_0 = 0.5;
    HMat = [NormalH_0, CompensationH_0];
    theta_0_Mat = atan(HMat);
    
    theta_m_Mat = atan(sqrt((1 + HMat.^2).^(1/3) - 1));
    U_snap_Mat = HMat - sqrt((1 + HMat.^2).^(1/3) - 1);
    F_snap_Mat = 2 * (sin(theta_m_Mat) - tan(theta_m_Mat) .* cos(theta_0_Mat));
    

    CorrK = F_snap_Mat(1) / F_snap_Mat(2);
    RealE = (CorrK - 1) * CompStatus .* NormE + NormE;

    [~, SortIndex] = sort(NormE);
    UnitSum = size(NormE, 2);
    HeapPos = sym(zeros(1, UnitSum));

    for i = 1: UnitSum
        IfThisComp = CompStatus(SortIndex(i));
        U_all = sym(zeros(1, UnitSum));
        for j = 1: UnitSum
            if (j == i)
                U_all(j) = U_snap_Mat(IfThisComp + 1);
                continue;
            end
            IfThatComp = CompStatus(SortIndex(j));
            IsDown = (j < i);
            if (j == i + 1)
                U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), 0, -2);
            else
                U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), IsDown, 2);
            end
        end
        HeapPos(i) = sum(U_all);
    end
end