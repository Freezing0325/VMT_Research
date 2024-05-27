function [RealE, HeapPos] = VMT_CalHeapPos(NormE, CompStatus, CalMethod, ActiveStatus)
%VMT_CalHeapPos     计算当前一侧配置的串联VMT所对应的峰值位置。
%
%   NormE           所有单元的归一化刚度
%   CompStatus      补偿情况，-1表示增加性补偿，0表示无补偿，1表示减缩性补偿
%   CalMethod       计算算法，1：求解非线性方程；2：线性模型：泰勒级数2阶拟合；3：非线性模型：泰勒级数3阶拟合。
%   ActiveStatus    单元激活情况，只有一行，不输入时，认为所有单元均激活。

    NormalH_0 = 0.5;
    CompMinusH_0 = 0.375;
    CompPlusH_0 = 0.625;
    
    HMat = [CompPlusH_0, NormalH_0, CompMinusH_0];

    [F_snap_Mat, U_snap_Mat] = VMT_SingleGetFm(1, HMat, CalMethod);

    CorrE = F_snap_Mat(2) ./ F_snap_Mat;
    RealE = NormE .* CorrE(CompStatus + 2);
    
    UnitSum = size(NormE, 2);
    TempHeapPos = zeros(1, UnitSum);
    if (~exist('ActiveStatus', 'var') || size(ActiveStatus, 1) == 0)
        ActiveStatus = ones(1, UnitSum);
    end
    for i = 1: UnitSum
        if (~ActiveStatus(i))
            continue;
        end
        ThisComp = CompStatus(i);
        U_all = zeros(1, UnitSum);
        for j = 1: UnitSum
            if (j == i)
                U_all(j) = U_snap_Mat(ThisComp + 2);
                continue;
            end
            ThatComp = CompStatus(j);
            IsDown = (j < i) || (~ActiveStatus(j));
            U_all(j) = VMT_SingleGetU(RealE(j), HMat(ThatComp + 2), RealE(i) * F_snap_Mat(ThisComp + 2), IsDown, CalMethod);
            if (~isreal(U_all(j)))
                U_all(j) = VMT_SingleGetU(RealE(j), HMat(ThatComp + 2), RealE(i) * F_snap_Mat(ThisComp + 2), IsDown, -CalMethod);
            end
        end
        TempHeapPos(i) = sum(U_all);
    end
    HeapPos = TempHeapPos(TempHeapPos > 0);
end