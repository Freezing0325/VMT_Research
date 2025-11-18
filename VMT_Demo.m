%% 可重构性的验证
% LeftNormE = [1.0000  1.1424  2.8561  3.0019];
% RightNormE = [1.0000  2.1428  2.2485  5.9393];
LeftNormE = [1.0000  1.4528  4.3279  4.4279  6.3568];
RightNormE = [1.0000  1.8311  1.9311  4.9533  5.0533];
LeftComp = [0 0 1 0 1];
RightComp = [0 1 0 1 0];
OriginStatus = 0;
CalMethod = 2;
[AllSequence, BestInactive_Sequence] = VMT_GetAllPossibleSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, 3, []);
% 0111最佳：4514
% [AllSequence2, BestInactive_Sequence2] = VMT_GetAllPossibleSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, 3, [1 1 1]);
% 0100最佳：2415
% [AllSequence2, BestInactive_Sequence2] = VMT_GetAllPossibleSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, 3, [1 0 0]);
% 0110最佳：3414
% [AllSequence, BestInactive_Sequence] = VMT_GetAllPossibleSequence(LeftNormE, RightNormE, LeftComp, RightComp, OriginStatus, CalMethod, 3, [1 1 0]);