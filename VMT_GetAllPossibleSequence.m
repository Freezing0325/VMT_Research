function [AllSequence, BestInactive_Sequence] = VMT_GetAllPossibleSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, ActiveSum, FindSequence)
%VMT_GetAllPossibleSequence         获得已有设计通过预配置可以获得的所有可能序列
%
%   LeftNormE       左侧单元的归一化刚度
%   RightNormE      右侧单元的归一化刚度
%   LeftComp        左侧补偿情况
%   RightComp       右侧补偿情况
%   OriginStatus    初始状态
%   CalMethod       计算算法，1：求解非线性方程；2：线性模型：泰勒级数2阶拟合；3：非线性模型：泰勒级数3阶拟合。
%   ActiveSum       两侧单元激活数
%   FindSequence    待查找的序列



    UnitSum = size(LeftNormE, 2);
    InactiveSum = UnitSum - ActiveSum;

    AllSequence = zeros(1, 2^(ActiveSum-1));
    BestInactive_Sequence = zeros(2^(ActiveSum-1), 2 * InactiveSum + 1);
    BestInactive_Sequence(:, 1) = 1;
    ChangeRange = ones(InactiveSum * 2, 1) * [2, 1, UnitSum];
    AllDesignSum = 0;
    OKDesignSum = 0;
    if (size(FindSequence, 1) == 0 || ~exist('FindSequence', 'var'))
        FindSequence = zeros(1, ActiveSum); % [1 0 1 0 1 0 1];
    end
    
    % NumOfVariables = size(ChangeRange, 1);
    [~, ChangeSum] = GetCombination(ChangeRange, 0, []);
    ThisCombination = [];
    for ChangeNo = 1: ChangeSum
        
        [ThisCombination, ~] = GetCombination(ChangeRange, ChangeNo, ThisCombination);
        
        % 要求前InactiveSum个数必须是升序的，后InactiveSum个数也是升序的。
        SortFlag = true;
        for i = 1: InactiveSum - 1
            if (ThisCombination(i) >= ThisCombination(i + 1) - 10e-6 || ...
                    ThisCombination(i + InactiveSum) >= ThisCombination(i + InactiveSum + 1) - 10e-6)
                SortFlag = false;
                break;
            end
        end
        if (~SortFlag)
            continue;
        end

        AllDesignSum = AllDesignSum + 1;
        ThisActiveStatus = ones(2, UnitSum);
        ThisActiveStatus(1, ThisCombination(1: InactiveSum)) = 0;
        ThisActiveStatus(2, ThisCombination(InactiveSum + 1: 2 * InactiveSum)) = 0;
        [ThisSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, ThisActiveStatus);
        if (MaxForceDiff < 1 && ThisSequence(1) == OriginStatus)
            OKDesignSum = OKDesignSum + 1;
            if (all(FindSequence == ThisSequence))
                for i = 1: InactiveSum * 2
                    fprintf('%d ', ThisCombination(i));
                end
                fprintf('\n');
            end
            ThisSequenceNo = 0;
            for i = 1: ActiveSum
                ThisSequenceNo = ThisSequenceNo + ThisSequence(i) * 2^(ActiveSum - i);
            end
            RealNo = ThisSequenceNo -  2^(ActiveSum-1) + 1;
            AllSequence(RealNo) = AllSequence(RealNo) + 1;
            if (MaxForceDiff < BestInactive_Sequence(RealNo, 1))
                BestInactive_Sequence(RealNo, 1) = MaxForceDiff;
                BestInactive_Sequence(RealNo, 2: 2 * InactiveSum + 1) = ThisCombination';
            end
        end
    end
end