OriginStatus = 1;   

LeftNormE = [1  4.5  5];
RightNormE = [1  3.5  6];

GoalSequence = [1 0 1];
StepNum = size(GoalSequence, 2);
ChangeInfo = GoalSequence - [OriginStatus, GoalSequence(1: StepNum - 1)];
LeftComp = -ChangeInfo;
RightComp = ChangeInfo;

[RealE_L, HeapPos_L] = VMT_CalHeapPos(LeftNormE, LeftComp, 3, []);
[RealE_R, HeapPos_R] = VMT_CalHeapPos(RightNormE, RightComp, 3, []);

