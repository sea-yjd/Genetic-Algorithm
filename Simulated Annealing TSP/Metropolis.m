function [S,R] = Metropolis(S1,S2,D,T)
%% Metropolis准则函数
% input:
% S1   当前解
% S2   新解
% D    距离矩阵（两两城市间距离）
% T    当前温度
% output：
% S    下一个当前解
% R    下一个当前解的路线距离

R1 = PathLength(D,S1);   %计算路线长度
N = length(S1);   % 得到城市的个数

R2 = PathLength(D,S2);   %计算路线长度
dC = R2 - R1;   %计算能量之差
if dC < 0   %如果能量降低，接受新路线
    S = S2;
    R = R2;
elseif exp(-dC/T) >= rand    % 以exp(-dC/T)概率接受新路线
    S = S2;
    R = R2;
else
    S = S1;
    R = R1;
end
