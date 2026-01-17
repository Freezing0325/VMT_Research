function [LeftComp, RightComp] = VMT_GetComp(GoalSequence, OriginStatus)
    GoalSequence_Hat = [OriginStatus, GoalSequence(1: end-1)];
    CompSide = GoalSequence - GoalSequence_Hat; 
    LeftComp = CompSide == -1;
    RightComp = CompSide == 1;
end
