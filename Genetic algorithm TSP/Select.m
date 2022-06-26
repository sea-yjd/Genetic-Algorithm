function SelCh = Select(Chrom,FitnV,GGAP)
%% 选择操作
% 输入：
% Chrom 种群
% FitnV 适应度值
% GGAP  选择概率
% 输出：
% SelCh  被选择的个体

NIND = size(Chrom,1);
NSel = max(floor(NIND*GGAP + .5), 2);
ChrIx = Sus(FitnV,NSel);
SelCh = Chrom(ChrIx,:);
