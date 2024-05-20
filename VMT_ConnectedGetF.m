function F = VMT_ConnectedGetF(E, H_0, U, Status, CalMethod)
%VMT_ConnectedGetU   对于串联的多个VMT单元，根据其总位移及双稳态状态计算受力
%
%   E               刚度矩阵
%   H_0             初始高度矩阵
%   U               总位移
%   Status          双稳态状态矩阵，0为上凸，1为下凸
%   CalMethod       计算方法，1：非线性精确求解，
%                   2：线性小变形单元（弹簧），二阶泰勒级数拟合
%                   3：线弹性大变形单元（均匀），三阶泰勒级数拟合
UnitSum = size(E, 2);
U_all = sym('U_', [1, UnitSum]);
syms F_all;

switch CalMethod
    case 2
        for i = 1: UnitSum
            U_all(i) = VMT_SingleGetU(E(i), H_0(i), F_all, Status(i), CalMethod);
        end
        eqn = sum(U_all) == U;
        F = solve(eqn, F_all, 0);
end
end