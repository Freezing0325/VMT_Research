function [Fm, Um] = VMT_SingleGetFm(EA, H_0, CalMethod)
%VMT_SingleGetFm   对于单个VMT单元，根据其刚度与初始高度计算极大值力
%
%   EA              抗拉刚度
%   H_0             初始高度
%   CalMethod       计算方式
%
%   输出：
%   Fm              极大值力
%   Um              极大值力对应的位移

    switch CalMethod
        case {1, 3}
            Fm = EA .* 2/(3*sqrt(3)) .* H_0.^3 ./ sqrt(H_0.^2 + 1);
            Um = H_0 - H_0 ./ sqrt(2 * H_0.^2 + 3);
        case 2
            Fm = 2 * EA .* (1 - (H_0.^2 + 1) .^ (-1/3)).^(3/2);
            Um = H_0 - ((H_0.^2 + 1) .^ (1/3) - 1) .^ (1/2);
    end
end