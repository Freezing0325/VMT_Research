function [RealE, HeapPos] = VMT_CalHeapPos_3(NormK, CompStatus)
%VMT_CalHeapPos_3     用三阶泰勒级数计算当前一侧配置的串联VMT所对应的峰值位置，用来符号求解的时候专用
%
%   NormK           所有单元的归一化刚度
%   CompStatus      补偿情况，1表示有补偿，0表示无补偿

    CompensationH_0 = 0.375;
    NormalH_0 = 0.5;
    HMat = [NormalH_0, CompensationH_0];
    
    U_snap_Mat = HMat - HMat ./ sqrt(2*HMat.^2 + 3);
    F_snap_Mat = 2 * HMat.^3 ./ (3 * sqrt(3) .* sqrt(HMat.^2 + 1));
    

    CorrK = F_snap_Mat(1) / F_snap_Mat(2);
    RealE = (CorrK - 1) * CompStatus .* NormK + NormK;

    [~, SortIndex] = sort(NormK);
    UnitSum = size(NormK, 2);
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
                U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), 0, -3);
            else
                U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), IsDown, 3);
            end
        end
        HeapPos(i) = sum(U_all);
    end
end