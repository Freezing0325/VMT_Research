function U = VMT_ConnectedGetU(E, H_0, F, Status, CalMethod)
%VMT_ConnectedGetU   对于串联的多个VMT单元，根据其受力大小及双稳态状态计算总位移
%
%   E               刚度矩阵
%   H_0             初始高度矩阵
%   F               受力
%   Status          双稳态状态矩阵，0为上凸，1为下凸
%   CalMethod       计算方法，1：非线性精确求解，
%                   2：线性小变形单元（弹簧），二阶泰勒级数拟合
%                   3：线弹性大变形单元（均匀），三阶泰勒级数拟合
U = 0;
UnitSum = size(E, 2);
for i = 1: UnitSum
    U = U + VMT_SingleGetU(E(i), H_0(i), F, Status(i), CalMethod);
end
end