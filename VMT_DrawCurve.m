function [U_Load, AllF] = VMT_DrawCurve(E, H_0, CalMethod)
%VMT_DrawCurve   对于串联VMT单元，根据其刚度和初始高度绘制力-位移曲线
%
%   E               刚度
%   H_0             初始高度
%   CalMethod       计算方法
%
%   输出：
%   U_Load          曲线的位移
%   AllF            曲线的力

    A = 1;    
    syms U_temp H_0_temp EA_temp;
    theta = atan(H_0_temp - U_temp);
    switch CalMethod
        case 1
            F = EA_temp * ((1 + H_0_temp^2) / (1 + (H_0_temp - U_temp)^2) - 1) * sin(theta);
        case 2
            F = 2 * EA_temp * (H_0_temp - U_temp) * (((H_0_temp - U_temp)^2 + 1)^(-1/2)-(H_0_temp^2+1)^(-1/2));
    end
    [Fm, MaxPt_all] = VMT_SingleGetFm(E * A, H_0, CalMethod);
    [~, ChangeIndex] = sort(Fm);
    
    UnitNum = size(E, 2);
    F_all = sym('F', [1, UnitNum]);
    U_all = sym('U', [1, UnitNum]);
    for i = 1: UnitNum
        F_all(i) = subs(F, {'U_temp', 'H_0_temp', 'EA_temp'}, {U_all(i), H_0(i), E(i) * A});
    end
    
    timestep = 0: 0.005: 1;
    U_Load = (sum(H_0) * 2 + 1) * timestep;
    CalStepSum = size(timestep, 2);
    
    
    AllF = zeros(CalStepSum, 1);
    DownNow = false;
    ChangeNum = 0;
    RandomSolveTimes = 0;
    
    uLocation = zeros(UnitNum, 2);
    uLocation(:, 1) = -0.25 * H_0';
    uLocation(:, 2) = MaxPt_all';
    
    sov_old = zeros(UnitNum, 1);
    eqn_all = sym('eqn', [1, UnitNum]);
    oldF = 0;
    
    for i = 1: UnitNum - 1
        eqn_all(i) = F_all(i) == F_all(i + 1);
    end

    for step = 1: CalStepSum
        ThisDelta = U_Load(step);
        eqn_all(UnitNum) = sum(U_all) == ThisDelta;
        sov = struct2cell(vpasolve(eqn_all, U_all, sov_old));
        sov = vertcat(sov{:});
        IsSolved = size(sov, 1);
        if (~IsSolved)
            if (~DownNow)
                ThisChangeIndex = ChangeIndex(ChangeNum + 1);
                if (ThisChangeIndex == 1)
                    uLocation(ThisChangeIndex, :) = [MaxPt_all(ThisChangeIndex), Inf];
                else
                    uLocation(ThisChangeIndex, :) = [2 * H_0(ThisChangeIndex) - MaxPt_all(ThisChangeIndex), Inf];
                end
            end
            RandomSolveTimes = RandomSolveTimes + 1;
            sov = struct2cell(vpasolve(eqn_all, U_all, uLocation));
            sov = vertcat(sov{:});
        end
        
        AllF(step) = double(subs(F_all(1), U_all(1), sov(1)));
        if (ChangeNum < UnitNum)
            ThisChangeIndex = ChangeIndex(ChangeNum + 1);
            if (~DownNow && AllF(step) < oldF)
                if (ThisChangeIndex == 1)
                    uLocation(ThisChangeIndex, :) = [MaxPt_all(ThisChangeIndex), Inf];
                else
                    uLocation(ThisChangeIndex, :) = [2 * H_0(ThisChangeIndex) - MaxPt_all(ThisChangeIndex), Inf];
                end
                ChangeNum = ChangeNum + 1;
                DownNow = true;
            end
            if (DownNow && AllF(step) > oldF)
                DownNow = false;
            end
        end
        oldF = AllF(step);
        sov_old = sov;
    end
    
    figure(1);
    plot(U_Load, AllF);
end
