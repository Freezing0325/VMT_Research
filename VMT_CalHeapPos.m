function [RealK, HeapPos] = VMT_CalHeapPos(NormK, CompStatus, CalMethod, ActiveStatus)
%VMT_CalHeapPos     计算当前一侧配置的串联VMT所对应的峰值位置。
%
%   NormK           所有单元的归一化刚度
%   CompStatus      补偿情况，1表示有补偿，0表示无补偿
%   CalMethod       计算算法，1：求解非线性方程；2：线性模型：泰勒级数2阶拟合；3：非线性模型：泰勒级数3阶拟合。
%   ActiveStatus    单元激活情况，只有一行，不输入时，认为所有单元均激活。
    CompensationH_0 = 0.375;
    NormalH_0 = 0.5;
    HMat = [NormalH_0, CompensationH_0];
    HMax_Mat = HMat ./ sqrt(2 * HMat.^2 + 3);
    theta_0_Mat = atan(HMat);
    if (CalMethod == 1)
        U_snap_Mat = HMat - HMat ./ sqrt(2*HMat.^2 + 3);
        F_snap_Mat = 2 * HMat.^3 ./ (3 * sqrt(3) .* sqrt(HMat.^2 + 1));
        syms H H_0;
        F = H / sqrt(H^2 + 1) * ((H_0^2 + 1) / (H^2 + 1) - 1);
    elseif (CalMethod == 2)
        theta_m_Mat = atan(sqrt((1 + HMat.^2).^(1/3) - 1));
        U_snap_Mat = HMat - sqrt((1 + HMat.^2).^(1/3) - 1);
        F_snap_Mat = sin(theta_m_Mat) - tan(theta_m_Mat) .* cos(theta_0_Mat);
    elseif (CalMethod == 3)
        U_snap_Mat = HMat - HMat ./ sqrt(2*HMat.^2 + 3);
        F_snap_Mat = 2 * HMat.^3 ./ (3 * sqrt(3) .* sqrt(HMat.^2 + 1));
    end

    CorrK = F_snap_Mat(1) / F_snap_Mat(2);
    RealK = (CorrK - 1) * CompStatus .* NormK + NormK;

    [~, SortIndex] = sort(NormK);
    UnitSum = size(NormK, 2);
    TempHeapPos = zeros(1, UnitSum);
    if (~exist('ActiveStatus', 'var') || size(ActiveStatus, 1) == 0)
        ActiveStatus = ones(1, UnitSum);
    end
    for i = 1: UnitSum
        if (~ActiveStatus(i))
            continue;
        end
        IfThisComp = CompStatus(SortIndex(i));
        U_all = zeros(1, UnitSum);
        for j = 1: UnitSum
            if (j == i)
                U_all(j) = U_snap_Mat(IfThisComp + 1);
                continue;
            end
            IfThatComp = CompStatus(SortIndex(j));
            IsDown = (j < i) || (~ActiveStatus(j));
            if (CalMethod == 1)
                ThisF = subs(F, H_0, HMat(IfThatComp + 1));
                SolutionPos = [-Inf, -HMax_Mat(IfThatComp + 1)] * (IsDown * 2 - 1);
                U_all(j) = HMat(IfThatComp + 1) - vpasolve(ThisF == RealK(SortIndex(i)) / RealK(SortIndex(j)) * F_snap_Mat(IfThisComp + 1), SolutionPos);
            elseif (CalMethod == 2)
                U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), IsDown, 2);
                if (~isreal(U_all(j)))
                    U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), IsDown, -2);
                end
            elseif (CalMethod == 3)
                U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), IsDown, 3);
                if (~isreal(U_all(j)))
                    U_all(j) = VMT_SingleGetU(RealE(SortIndex(j)), HMat(IfThatComp + 1), RealE(SortIndex(i)) * F_snap_Mat(IfThisComp + 1), IsDown, -3);
                end
            end
        end
        TempHeapPos(i) = sum(U_all);
    end
    HeapPos = TempHeapPos(TempHeapPos > 0);
end