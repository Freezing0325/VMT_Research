

OriginStatus = 1;   

LeftNormK = [1  1.7  1.8  2.6  2.7  3.6  3.7  4.7  4.8];
RightNormK = [1  1.5  2  2.5  3  3.5  4  4.5  5];

GoalSequence = [1 0 1 0 1 0 1 0 1];   % 1：左高右低，0：左低右高。
StepSum = size(GoalSequence, 2);    % 序列长度
CompSide = GoalSequence - [OriginStatus, GoalSequence(1: StepSum - 1)]; % 需要施加补偿的一侧，0：不需要，-1：左侧，1：右侧。
ActiveStatus = [1 1 1 1 1 1 1 1 1  ; ...
                1 1 1 1 1 1 1 1 1];

ActiveSum = 7;
[AllSequence, BestInactive_Sequence] = VMT_GetAllPossibleSequence(LeftNormK, RightNormK, CompSide, OriginStatus, ActiveSum, []);
DiffSequenceSum = size(find(AllSequence), 2);

fprintf('预压缩单元数量：%d，实现的序列种数：%d\n', StepSum - ActiveSum, DiffSequenceSum);

%%

RightNormE = [1 1.5 2 2.5 3 3.5 4 4.5 5];
GoalSequence = [1 0 1 0 1 0 1 0 1];   % 1：左高右低，0：左低右高。
OriginStatus = 1;
BeginE = [];
% BeginE = [1  1.9  2  2.9  3  4.4  4.5  6.9  7];
[BestRealE_Left, BestE, RealE_Right] = VMT_DesignOtherSide(RightNormE, GoalSequence, OriginStatus, BeginE);


