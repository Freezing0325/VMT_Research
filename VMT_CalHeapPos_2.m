function [RealE, HeapPos] = VMT_CalHeapPos_2(NormE, CompStatus)
%VMT_CalHeapPos_2     用二、三阶泰勒级数、线性模型计算当前一侧配置的串联VMT所对应的峰值位置
%
%   NormK           所有单元的归一化刚度
%   CompStatus      补偿情况，1表示有补偿，0表示无补偿

    NormalH_0 = 0.5;
    CompMinusH_0 = 0.375;
    CompPlusH_0 = 0.625;
    
    HMat = [CompPlusH_0, NormalH_0, CompMinusH_0];
    
    [F_snap_Mat, U_snap_Mat] = VMT_SingleGetFm(1, HMat, 2);
    

    CorrE = F_snap_Mat(2) ./ F_snap_Mat;
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
            if (j == i + 1)
                U_all(j) = VMT_SingleGetU(RealE(j), HMat(IfThatComp + 2), RealE(i) * F_snap_Mat(IfThisComp + 2), 0, -2);
            else
                U_all(j) = VMT_SingleGetU(RealE(j), HMat(IfThatComp + 2), RealE(i) * F_snap_Mat(IfThisComp + 2), IsDown, 2);
            end
        end
        HeapPos(i) = sum(U_all);
    end
end