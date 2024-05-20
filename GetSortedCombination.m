function AllCombination = GetSortedCombination(VaryRange, IfRigid)
%GetCombination     获得顺序排列的所有组合
%   VaryRange存储各个变量的变化范围，
%   IfRigid表示是否严格大于。
%   其中每一行代表一个变量的变化，第一个元素为起始位置，第二个元素为步长，第三个元素为结束位置，变化从起始位置开始，若大于结束位置则结束。
%   注意，所有组合的排序方式为：从后面的变量开始逐个变化，遍历完一个变量后向前"进位"。
%   输出所有的组合方式，每一行为一种组合方式

    NumOfVariables = size(VaryRange, 1);

    AllCombination = [];
    if (NumOfVariables == 1)
        AllCombination = (VaryRange(1, 1): VaryRange(1, 2): VaryRange(1, 3))';
        return;
    end
    for FirstNum = VaryRange(1, 1): VaryRange(1, 2): VaryRange(1, 3)
        if (IfRigid)
            NextBeginNum = VaryRange(2, 1) + max(floor((FirstNum - VaryRange(2, 1)) / VaryRange(2, 2) + 10e-6) + 1, 0) * VaryRange(2, 2);
        else
            NextBeginNum = VaryRange(2, 1) + max(ceil((FirstNum - VaryRange(2, 1)) / VaryRange(2, 2) - 10e-6), 0) * VaryRange(2, 2);
        end
        
        NextVaryRange = VaryRange(2: NumOfVariables, :);
        NextVaryRange(1, 1) = NextBeginNum;
        NextAllCombination = GetSortedCombination(NextVaryRange, IfRigid);
        NextVarySum = size(NextAllCombination, 1);
        AllCombination = [AllCombination; ones(NextVarySum, 1) * FirstNum, NextAllCombination];
    end
end


