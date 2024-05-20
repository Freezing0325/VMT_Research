function U = VMT_SingleGetU(E, H_0, F, Status, CalMethod)
%VMT_SingleGetU   对于单个VMT单元，根据其受力大小及双稳态状态计算其位移
%
%   E               刚度
%   H_0             初始高度
%   F               受力
%   Status          双稳态状态，0为上凸，1为下凸
%   CalMethod       计算方法，1：非线性精确求解，
%                   2：线性小变形单元（弹簧），二阶泰勒级数拟合，-2：线性小变形单元（弹簧），力F接近极值，三阶泰勒级数拟合。
%                   3：线弹性大变形单元（均匀），三阶泰勒级数拟合，-3：线弹性大变形单元（均匀），力F接近极值

switch CalMethod
    case 1
        H_m = H_0 / sqrt(2 * H_0^2 + 3);
        syms H;
        if (Status == 0)
            solve_pos = [H_m, inf];
        else
            solve_pos = [-inf, -H_m];
        end
        H_solve = vpasolve(F == E * H / (sqrt(H^2 + 1)) * ((H_0^2 + 1) / (H^2 + 1) - 1), H, solve_pos);
        U = H_0 - H_solve;
    case -2
        theta_m = atan(sqrt((1 + H_0^2)^(1/3) - 1));
        Fm = 2 * (sin(theta_m) - tan(theta_m) .* cos(atan(H_0)));
        Eqn_a = 2 * (-(2/(H_0^2 + 1)^(5/6) - 5/(2 * (H_0^2 + 1)^(7/6))));
        Eqn_b = 2 * 3 * ((H_0^2 + 1)^(1/3) - 1)^(1/2) / (2 * (H_0^2 + 1)^(5/6));
        Eqn_c = 0;
        Eqn_d = F / E - Fm; %RealE(SortIndex(i)) / RealE(SortIndex(j)) * F_snap_Mat(IfThisComp + 1) - Fm;
        Eqn_p = Eqn_c / Eqn_a - (Eqn_b / Eqn_a)^2 / 3;
        Eqn_q = Eqn_d / Eqn_a + 2 * (Eqn_b / (3 * Eqn_a))^3  - Eqn_b * Eqn_c / (3 * Eqn_a^2);
        Eqn_r = sqrt(-(Eqn_p/3)^3);
        Eqn_theta = acos(-Eqn_q / (2 * Eqn_r)) / 3;
        U = H_0 - sqrt((1 + H_0^2).^(1/3) - 1) -(2 * Eqn_r^(1/3) * cos(Eqn_theta) - Eqn_b / (3 * Eqn_a));
    case 2
        Eqn_a = 2 * (-(-1.5 * H_0 ./ (H_0.^2 + 1) .^ (5/2)));
        Eqn_b = 2 * (-(-H_0.^2 ./ (H_0.^2 + 1) .^ (3/2)));
        Eqn_c = -F / E * (Status * 2 - 1);
        Eqn_Delta = Eqn_b^2 - 4 * Eqn_a * Eqn_c;
        U = -(sqrt(Eqn_Delta) - Eqn_b) / (2 * Eqn_a);
        if (Status == 1)
            U = 2 * H_0 - U;
        end

    case -3
        Eqn_d = -F / E * (Status * 2 - 1);
        H_m = H_0 / sqrt(2 * H_0^2 + 3);
        Fm = 2 * H_0^3 / (3 * sqrt(3) * sqrt(H_0^2 + 1));
        U = H_0 - H_m - sqrt((Fm - Eqn_d) * (3 * H_0^2 + 3)^(5/2) / (H_0 * (2 * H_0^2 + 3)^3));
    case 3
        Eqn_a = -(4*H_0^4 - 10*H_0^2 + 1) / (H_0^2 + 1)^(7/2);
        Eqn_b = -(3*H_0*(H_0^2 - 1)) / (H_0^2 + 1)^(5/2);
        Eqn_c = -2*H_0^2 / (H_0^2 + 1)^(3/2);
        Eqn_d = -F / E * (Status * 2 - 1);
        Eqn_p = Eqn_c / Eqn_a - (Eqn_b / Eqn_a)^2 / 3;
        Eqn_q = Eqn_d / Eqn_a + 2 * (Eqn_b / (3 * Eqn_a))^3  - Eqn_b * Eqn_c / (3 * Eqn_a^2);
        Eqn_r = sqrt(-(Eqn_p/3)^3);
        Eqn_theta = acos(-Eqn_q / (2 * Eqn_r)) / 3;
        U = 2 * Eqn_r^(1/3) * cos(Eqn_theta - 2 * pi / 3) - Eqn_b / (3 * Eqn_a);
        if (Status == 1)
            U = 2 * H_0 - U;
        end
end


end