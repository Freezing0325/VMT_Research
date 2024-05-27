OriginStatus = 1;   

LeftNormE = [1  3  4];
RightNormE = [1  2  5];

GoalSequence = [1 0 1];
StepNum = size(GoalSequence, 2);
ChangeInfo = GoalSequence - [OriginStatus, GoalSequence(1: StepNum - 1)];
LeftComp = -ChangeInfo;
RightComp = ChangeInfo;

[RealE_L, HeapPos_L] = VMT_CalHeapPos(LeftNormE, LeftComp, 3, []);
[RealE_R, HeapPos_R] = VMT_CalHeapPos(RightNormE, RightComp, 3, []);

[PredSequence, MaxForceDiff] = VMT_GetSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, 3, []);