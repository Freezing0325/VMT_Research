function [ThisCombination, VarySum] = GetCombination(VaryRange, No, LastCombination)
%GetCombination     获得某个组合
%   VaryRange存储各个变量的变化范围，
%   其中每一行代表一个变量的变化，第一个元素为起始位置，第二个元素为步长，第三个元素为结束位置，起始位置和结束位置都会取到。
%   No表示想要取出的组合序号。
%   注意，所有组合的排序方式为：从后面的变量开始逐个变化，遍历完一个变量后向前"进位"。

    GetSum = No == 0;
    ExistLastCombination = size(LastCombination, 1) ~= 0;
    VarySum = 0;
    NumOfVariables = size(VaryRange, 1);

    if (ExistLastCombination)
        ThisCombination = LastCombination;
        
        for i = NumOfVariables: -1: 1
            ThisCombination(i) = ThisCombination(i) + VaryRange(i, 2);
            if (ThisCombination(i) > VaryRange(i, 3) + 10e-6)
                ThisCombination(i) = VaryRange(i, 1);
                continue;
            end
            break;
        end
    else
        
        VaryNumOfEveryVariable = floor((VaryRange(:, 3) - VaryRange(:, 1)) ./ VaryRange(:, 2)) + 1;
    
        if (GetSum)
            VaryNo = ones(NumOfVariables + 1, 1);
            for i = NumOfVariables: -1: 1
                VaryNo(i) = VaryNo(i + 1) * VaryNumOfEveryVariable(i);
            end
            VarySum = VaryNo(1);
            ThisCombination = [];
        else
            
            ThisCombination = VaryRange(:, 1);
            ThisNo = No;
            for i = NumOfVariables: -1: 1
                ThisCombination(i) = ThisCombination(i) + mod(ThisNo - 1, VaryNumOfEveryVariable(i)) * VaryRange(i, 2);
                ThisNo = floor((ThisNo - 1) / VaryNumOfEveryVariable(i));
                if (ThisNo == 0)
                    break;
                end
            end
        end
    end
end


