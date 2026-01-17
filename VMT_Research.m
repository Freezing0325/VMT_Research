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
% 所有两行，上一行为初始，下一行为最终

% B
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.661	1.905	2.334	6.755	10.823 12.826], [1.661	1.905	2.334	4.017	7.878 9.320], [0 1 1 1 1 1], [0 1 1 0 1 1], 1, 2, [])
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.661	1.905	2.334	5.02	7.878 9.320]/1.545, [1.547	1.905	2.334	6.755	10.823 12.826]  /1.545, [0 1 1 0 1 1], [0 1 1 1 1 1], 0, 2, []);
VMT_ReportConfig(1,[1	1.905	2.334	(2.490+6.441)/2	7.878 9.320 1	1.905	2.334	6.755	10.823 12.826]/1.547, [0 0 0 1 1 1; 0 1 1 0 1 1;0 1 1 1 1 1],0, 2);

% R
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1	1.477	4.691	6.937	9.07	9.798], [1	1.477	2.085	6.937	9.32	12.826], [0 1 1 1 1 0], [0 1 0 1 1 1], 1, 2, [])
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1	1.477	4.691	7.879	9.07	9.798]/1.545, [1	1.477	2.085	5.814	6.937	12.826]/1.545, [0 1 1 1 1 0], [0 1 0 1 1 1], 1, 2, [])
% VMT_ReportConfig(1,[1	1.477	4.691	6.937	9.07	9.798  1	1.477	3.093	5.814	6.937	12.826]/1.545, [1 1 0 0 0 1; 0 1 1 1 1 0;0 1 0 1 1 1],1 ,2);

% G
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1	1.074	3.033	3.7595	7.1595	8.293], [1.074	4.486	5.256	5.8645	6.3355	7.1595], [0 0 1 1 1 1], [0 1 1 1 0 1], 0, 2, []);
% VMT_ReportConfig(1, [1	2.598	3.033	3.7595	7.1595	8.293   (0.646+1.074)/2	4.486	5.256	5.8645	6.3355	7.1595], [0 1 1 1 0 0;0 0 1 1 1 1;0 1 1 1 0 1], 0, 2); % 实验很难做
% VMT_ReportConfig(1, [1	(2.598+1.348)/2	3.033	3.7595	7.1595	8.293   (0.646+1)/2	4.486	5.256	5.8645	6.3355	7.1595], [0 1 1 1 0 0;0 0 1 1 1 1;0 1 1 1 0 1], 0, 2);
% VMT_ReportConfig(1, [1	(2.598+1.305)/2	3.033	3.7595	7.1595	8.293   0.646	3.486	4.256	5.8645	6.3355	7.1595], [0 1 1 1 0 0;0 0 1 1 1 1;0 1 1 1 0 1], 0, 2);


% 10101010 Most Complex Rule
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.477	1.547	3.513	3.940	8.129	8.246	12.826], [1	2.334	2.49	5.814	6.023	8.759	9.121], [1 0 1 0 1 0 1], [0 1 0 1 0 1 0], 1, 2, [])


% Rule 27
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.000 	3.033 	3.894 	4.486 	6.026 	6.998 	8.293], [1.348	1.61	4.368	4.486	4.5525	5.703	8.293], [0 1 0 1 1 1 1], [0 0 1 1 0 1 1], 1, 2, [])

%[PredSequence, MaxForceDiff] = VMT_GetSequence([1.000 	2.422 	2.598 	3.033 	6.026 	6.998 	8.293], [1.348	1.61	4.368	4.486	4.553	5.703	7.160], [0 1 0 1 1 1 1], [0 0 1 1 0 1 1], 1, 2, [])
% [PredSequence, MaxForceDiff] = VMT_GetSequence([1.000 	2.422 	2.598 	3.033 	6.026 	6.998 	8.293], [1.348	1.61	4.368	4.486	4.924	5.703	7.160], [0 1 0 1 1 1 1], [0 0 1 1 0 1 1], 1, 2, [])
%[PredSequence, MaxForceDiff] = VMT_GetSequence([1.000 	2.422 	2.598 	3.033 	6.026 	6.998 	8.293], [1.6	1.61 4.368	4.486	4.924	5.703	7.160], [0 1 0 1 1 1 1], [0 0 1 1 0 1 1], 1, 2, [])
%

%% 4阶段全部变形序列
OriginStatus = 0;
CalMethod = 2;
% 0000 第一阶段有问题
% LeftNormE = [0.647, 1.348, 1.610, 3.484];
% RightNormE = [0.647, 1.074, 1.346, 3.078];
% GoalSequence = [0 0 0 0];

% 0001 第一阶段有问题
% LeftNormE = [0.646 1.348 1.610 3.078];
% RightNormE = [0.646 1.074 1.346 5.703];
% GoalSequence = [0 0 0 1];

% 0010 第一阶段有问题
% LeftNormE = [0.646 1.346 2.363 5.703];
% RightNormE = [0.646 1.074 3.033 3.484];
% GoalSequence = [0 0 1 0];

% 0011
% LeftNormE = [0.646 1.610 2.363 2.598];
% RightNormE = [0.646 1 3.033 3.894];
% GoalSequence = [0 0 1 1];

% 0100 最终位移不太好
% LeftNormE = [0.646 1 3.033 3.484];
% RightNormE = [0.646 1.509 2.597 3.078];
% GoalSequence = [0 1 0 0];

% 0110 最终位移不太好
% LeftNormE = [0.646 1 1.610 5.073];
% RightNormE = [0.646 1.509 2.336 2.597];
% GoalSequence = [0 1 1 0];

% 0111 第二段-3警告
% LeftNormE = [0.646 (1.346+2.336)/2 2.363 4.165];
% RightNormE = [0.646 3.033 3.484 6.769];
% GoalSequence = [0 1 1 1];
% [PredSequence, MaxForceDiff] = VMT_GetSequence([0.646 1 2.363 2.598], [0.646 3.033 3.078 6.769], [0 0 0 0], [0 1 0 0], OriginStatus, CalMethod, [])
% [PredSequence, MaxForceDiff] = VMT_GetSequence([0.646 1 2.363 2.598], [0.646 1.509 3.078 6.769], [0 0 0 0], [0 1 0 0], OriginStatus, CalMethod, [])

[LeftComp, RightComp] = VMT_GetComp(GoalSequence, OriginStatus);
[PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, [])
VMT_ReportConfig(1,[LeftNormE RightNormE], [GoalSequence; LeftComp;RightComp],OriginStatus, CalMethod);
%% 9阶段超材料
LeftNormE =[1   1.1314    3.6175    3.7243    4.6044    4.7107    5.7353    5.8426    7.3346];
RightNormE = [1  1.3906    1.4911    4.0308    4.1332    5.0549    5.1577    6.2540    6.3573];
LeftComp = [0 0 1 0 1 0 1 0 1];
RightComp = [0 1 0 1 0 1 0 1 0];
OriginStatus = 0;
CalMethod = 2;
[PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, [])
%% 可调逻辑门 & 多路解码器 （会有-3警告）
% 1→1→0
LeftNormE =[1   (1.509+3.033)/2 8.293];
RightNormE = [1  (1.850+3.033)/2 (2.363+3.484)/2];
LeftComp = [0 1 1];
RightComp = [0 1 0];
OriginStatus = 1;
CalMethod = 2;
% [PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, [])
VMT_ReportConfig(1,[LeftNormE RightNormE], [1 1 0; LeftComp;RightComp],OriginStatus, CalMethod);
% 0→1→1
LeftNormE =[0.647   (1.348+1.610)/2 (6.026+3.033)/2];
RightNormE = [0.647  6.026 8.293];
% LeftNormE =[0.647   (1.948+1.610)/2 3.033];
% RightNormE = [0.647  3.033 8.293];
LeftComp = [0 0 1];
RightComp = [0 1 1];
OriginStatus = 0;
% [PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, [])
VMT_ReportConfig(1,[LeftNormE RightNormE], [0 1 1; LeftComp;RightComp],OriginStatus, CalMethod);

%% FPGA实验
% 四阶段0→1→0→1 第四段-3警告
OriginStatus = 0;
CalMethod = 2;
LeftNormE = [(0.646+1.114)/2 (1.346+1)/2 3.033 (3.484+6.769)/2];
RightNormE = [1 (1.509+3.033)/2  2.336 (4.486+8.293)/2];
GoalSequence = [0 1 0 1];
[LeftComp, RightComp] = VMT_GetComp(GoalSequence, OriginStatus);
% [PredSequence, MaxForceDiff] = VMT_GetSequence([0.646 1 3.033 (3.484+6.769)/2], [1 (1.509+3.033)/2  2.336 (4.486+8.293)/2], [0 0 1 0], [0 1 0 1], 0, 2, [])
[PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod)
VMT_ReportConfig(1,[LeftNormE RightNormE], [GoalSequence; LeftComp;RightComp],OriginStatus, CalMethod);
%% FE FPGA
% 五阶段0→1→0→1→0
LeftNormE = [1 1.453 25.617/5.919 4.428 37.626/5.919];
RightNormE = [1 10.838/5.919  1.931 29.318/5.919 5.053];
LeftComp = [0 0 1 0 1];
RightComp = [0 1 0 1 0];
OriginStatus = 0;
CalMethod = 2;
ActiveSum = 3;
FindSequence = [];
[PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, []);
[AllSequence, BestInactive_Sequence] = VMT_GetAllPossibleSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, ActiveSum, FindSequence);
%% 上升段斜率的大致估计
VMT_Init();
global a Normal_h Output_h Output_a;
EA_ka = LeftNormE;
UnitSum = size(EA_ka, 2);
[RealEA_ka, HeapPos] = VMT_CalHeapPos(EA_ka, LeftComp, 2, []);
HeapPos = HeapPos;
Material_U0 = (0:UnitSum)*Normal_h/a*2;
Fm = VMT_SingleGetFm(EA_ka(1), Normal_h/a, 2);
EqStiffness = EA_ka/EA_ka(1)*Fm ./ (HeapPos - Material_U0(1:end-1));
F_out_true = 2 * VMT_SingleGetFm(2.8821 * Output_a, Output_h / Output_a, 2);
