function FitnV = Fitness(len)
%% 适应度函数
% 输入：
% len   个体的长度（TSP距离）
% 输出：
% FitnV  个体的适应度值

FitnV = 1./len;