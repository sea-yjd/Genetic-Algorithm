clc;
clear;
close all;
load citys_data.mat
tic
T0 = 1000;   %初始温度
Tend = 1e-3;   %终止温度
L = 200;     %各温度下的迭代次数（链长）
q = 0.985;     %降温速率
X = citys;
%%
D = Distanse(X);    % 计算距离矩阵
N = size(D,1);   %城市的个数
%% 计算迭代的次数Time
syms x;
eq = 1000*(q)^x == num2str(Tend);
Time = ceil(double(solve(eq,x)));
m_c = 50;
trace_new = zeros(Time,1);
%% 蒙特卡罗实验
for m = 1:m_c
    trace = [];
%% 初始解
S1 = randperm(N);      %随机产生一个初始路线

%% 画出随机解得路径图
% DrawPath(S1,X);
% pause(0.0001);
%% 输出随机解的路径和总距离
% disp('初始种群中的一个随机值：');
% OutputPath(S1);
% Rlength = PathLength(D,S1);
% disp(['总距离：',num2str(Rlength)]);

% Time = ceil(double(solve(['1000 * (0.9)^x =',num2str(Tend)])));
count = 0;           %迭代计数
Obj = zeros(Time,1);     %目标值矩阵初始化
track = zeros(Time,N);    %每代的最优路线矩阵初始化
T0 = 1000;   %初始温度
%% 迭代
while T0 > Tend
    count = count +1;    %更新迭代次数
    temp = zeros(L , N+1);
    for k = 1:L
        %% 产生新解
        S2 = NewAnswer(S1);
        %% Metropolis法则判断是否接受新解
        [S1,R] = Metropolis(S1,S2,D,T0);      %Metropolis抽样算法
        temp(k,:) = [S1 R];    %记录下一路线及其路程
    end
    %% 记录每次迭代过程的最优路线
    [d0,index] = min(temp(:,end));     %找出当前温度下最优路线
    if count==1 || d0 < Obj(count - 1)
        Obj(count) = d0;    %如果当前温度下最优路程小于上一路程，则记录当前路程
    else
        Obj(count) = Obj(count - 1);   %如果当前温度下最优路程大于上一路程，则记录上一路程
    end
    track(count,:) = temp(index,1:end-1);    %记录当前温度的最优路线
    T0 = q * T0;      %降温
%     fprintf(1,'%d\n',count);     %输出当前迭代次数

    trace = [trace;Obj(count)];
end
trace_new = trace_new + trace;
end
toc
%% 优化过程迭代图
figure 
plot(1:count,trace_new./m_c);
xlabel('代数');
ylabel('最优值/km');
axis([0,1000,14000,40000]);
title('模拟退火算法50次蒙特卡罗收敛曲线');
grid on;

%% 最优解的路径图
% DrawPath(track(end,:),X);
tourGbest = track(end,:);
% tourGbest = Chrom(minInd,:);
n = size(citys,1);
figure
hold on
plot([citys(tourGbest(1),1),citys(tourGbest(n),1)],[citys(tourGbest(1),2),...
    citys(tourGbest(n),2)],'ks-','Markersize',8,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
text(citys(1,1),citys(1,2),['   ' num2str(1)]);
hold on
for i=2:n
    plot([citys(tourGbest(i-1),1),citys(tourGbest(i),1)],[citys(tourGbest(i-1),2),...
        citys(tourGbest(i),2)],'ks-','Markersize',8,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
    text(citys(i,1),citys(i,2),['  ' num2str(i)]);
    hold on
end
legend('规划路径')
title([ '模拟退火算法优化路径（最短距离）：' num2str(PathLength(D,tourGbest))],'fontsize',10);
xlabel('城市位置横坐标/km','fontsize',10)
ylabel('城市位置纵坐标/km','fontsize',10)

grid on
x = citys(:,1);
y = citys(:,2);
startx=x(tourGbest(1)); %起点x坐标
starty=y(tourGbest(1)); %起点y坐标
endx=x(tourGbest(n));
endy=y(tourGbest(n));
text(startx,starty,'    起点'); %标记起点
text(endx,endy,'    终点')%标记终点
set(gca,'LineWidth',1.5);  %边框加粗,美观
axis([1000 1.1*max(x) 500 1.1*max(y)]); %设置尺寸大小

%% 输出最优解的路线和总距离
disp('最优解：');
S = track(end,:);
p = OutputPath(S);
disp(['总距离:',num2str(PathLength(D,S))]);
disp('============================');
















