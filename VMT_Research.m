close all;
U = -0.5: 0.01: 1.25;
H_0 = 0.375;
theta_0 = atan(H_0);
H = H_0 - U;
theta = atan(H);
EA_ka = 1;

%% 线性模型
% 在H_0 = 0.5和H_0 = 0.375的情况下，二阶拟合解的效果比三阶好。但在最大值处，三阶的效果比二阶好得多。
F = 2 * EA_ka * (sin(theta)/cos(theta_0) - tan(theta));
figure(1);
plot(U, F);
hold on;
U_app = -0.5: 0.01: 0.25;
Func_1 = -H_0^2/(H_0^2 + 1);
Func_2 = -3/2 * H_0 / (H_0^2 + 1)^2;
Func_3 = (4 * H_0^2 - 1) / (2 * (H_0^2 + 1)^3);
F_app_2 = 2 * EA_ka * (-Func_1 * U_app + Func_2 * U_app .^ 2);
F_app_3 = 2 * EA_ka * (-Func_1 * U_app + Func_2 * U_app .^ 2 - Func_3 * U_app .^ 3);

H_m = ((H_0^2 + 1)^(1/3) - 1)^(1/2);
F_m = ((H_0^2 + 1)^(1/3) - 1)^(3/2);
Func_m_2 = -3 * H_m/ (2 * (H_0^2 + 1)^(1/3));
Func_m_3 = 2 / (H_0^2 + 1)^(1/3) - 5 / (2 * (H_0^2 + 1)^(2/3));
U_m_app = -0.5: 0.01: 0.25;
H_m_app = H_0 - U_m_app - H_m;
F_m_app_2 = 2 * EA_ka * (F_m + Func_m_2 * H_m_app.^2);
F_m_app_3 = 2 * EA_ka * (F_m + Func_m_2 * H_m_app.^2 + Func_m_3 * H_m_app.^3);


plot(U_app, F_app_2);
plot(U_app, F_app_3);
plot(U_m_app, F_m_app_2);
plot(U_m_app, F_m_app_3);
legend('准确解', '二阶拟合解', '三阶拟合解', '最大值处二阶', '最大值处三阶', 'Location','southeast');
xlabel('U');
ylabel('F');
grid on;

% 对SingleGetU函数的验证
figure(2);
plot(U, F);
hold on;
PtNum = max(size(U_app));
MaxFPtNum = find(F(1:PtNum) == max(F(1: PtNum)));

U_Cal = zeros(1, MaxFPtNum);
for i = 1: MaxFPtNum
    U_Cal(i) = VMT_SingleGetU(EA_ka, H_0, F(i), 0, 2);
    if (~isreal(U_Cal(i)))
        U_Cal(i) = VMT_SingleGetU(EA_ka, H_0, F(i), 0, -2);
    end
end
plot(U_Cal, F(1: MaxFPtNum));
legend('准确解', '计算得解', 'Location','southeast');
xlabel('U');
ylabel('F');
grid on;

%% 非线性模型
F_nonlinear = EA_ka * sin(theta) .* ((cos(theta)).^2 / (cos(theta_0))^2 - 1);
figure(3);
plot(U, F_nonlinear);
hold on;
U_app = U(1): 0.01: 0.25;
Func_1_nonlinear = (2 * H_0^2)/(H_0^2 + 1)^(3/2);
Func_2_nonlinear = (3*H_0 * (H_0^2 - 1))/(H_0^2 + 1)^(5/2);
Func_3_nonlinear = (4*H_0^4 - 10 * H_0^2 + 1)/ (H_0^2 + 1)^(7/2);
Func_4_nonlinear = (10*H_0^5-45*H_0^3+15*H_0)/(2*(H_0^2+1)^(9/2));
F_app_nonlinear = EA_ka * (Func_1_nonlinear * U_app + Func_2_nonlinear * U_app .^ 2 + Func_3_nonlinear * U_app .^ 3 + Func_4_nonlinear * U_app .^ 4);
F_app_nonlinear_2 = EA_ka * (Func_1_nonlinear * U_app + Func_2_nonlinear * U_app .^ 2);
F_app_nonlinear_3 = EA_ka * (Func_1_nonlinear * U_app + Func_2_nonlinear * U_app .^ 2 + Func_3_nonlinear * U_app .^ 3);
plot(U_app, F_app_nonlinear);
plot(U_app, F_app_nonlinear_2);
plot(U_app, F_app_nonlinear_3);
% U_app2 = 0: 0.01: 1;
% F_app2 = 0.1180 * (0.5 - U_app2) - 0.5590 * (0.5 - U_app2) .^ 3;
% plot(U_app2, F_app2);
legend('准确解', '四阶拟合解', '二阶拟合解', '三阶拟合解', 'Location','southeast');
xlabel('U');
ylabel('F_{nonlinear}');
grid on;

% 对SingleGetU函数的验证
figure(4);
plot(U, F_nonlinear);
hold on;
PtNum = max(size(U_app));
MaxFPtNum = find(F_nonlinear(1:PtNum) == max(F_nonlinear(1: PtNum)));

U_Cal = zeros(1, MaxFPtNum);
for i = 1: MaxFPtNum
    U_Cal(i) = VMT_SingleGetU(EA_ka, H_0, F_nonlinear(i), 0, 3);
    if (~isreal(U_Cal(i)))
        U_Cal(i) = VMT_SingleGetU(EA_ka, H_0, F_nonlinear(i), 0, -3);
    end
end
plot(U_Cal, F_nonlinear(1: MaxFPtNum));
legend('准确解', '计算得解', 'Location','southeast');
xlabel('U');
ylabel('F');
grid on;

%% 
syms U;
Func_1_nonlinear = (2 * H_0^2)/(H_0^2 + 1)^(3/2);
Func_2_nonlinear = (3*H_0 * (H_0^2 - 1))/(H_0^2 + 1)^(5/2);
Func_3_nonlinear = (4*H_0^4 - 10 * H_0^2 + 1)/ (H_0^2 + 1)^(7/2);
Func_4_nonlinear = (10*H_0^5-45*H_0^3+15*H_0)/(2*(H_0^2+1)^(9/2));
F_app_sym = Func_1_nonlinear * U + Func_2_nonlinear * U^2 +  Func_3_nonlinear * U^3 + Func_4_nonlinear * U^4;
%%
% B
%[PredSequence, MaxForceDiff] = VMT_GetSequence([1.661	1.905	2.334	6.755	10.823 12.826], [1.661	1.905	2.334	4.017	7.878 9.320], [0 1 1 1 1 1], [0 1 1 0 1 1], 1, 2, [])
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.547	1.905	2.334	6.755	10.823 12.826], [1.661	1.905	2.334	5.02	8.820 9.320], [0 1 1 1 1 1], [0 1 1 0 1 1], 1, 2, [])

% R
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1	1.477	4.691	6.937	9.07	9.798], [1	1.477	2.085	6.937	9.32	12.826], [0 1 1 1 1 0], [0 1 0 1 1 1], 1, 2, [])
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1	1.477	4.691	7.879	9.07	9.798], [1	1.477	2.085	5.814	6.937	12.826], [0 1 1 1 1 0], [0 1 0 1 1 1], 1, 2, [])

% G
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.074	4.486	5.256	5.8645	6.3355	7.1595], [1	1.074	3.033	3.7595	7.1595	8.293], [0 1 1 1 0 1], [0 0 1 1 1 1], 1, 2, [])


