clear all
clear
clc
%% 蚁群算法解决旅行商问题（随机初始化每个蚂蚁最初所在城市，利用轮盘赌法（根据城市之间的转移概率）选择下一个访问城市，计算每个蚂蚁的路径长度
%% 更新信息素）
load citys_data.mat
tic
%% 计算城市间相互距离
n = size(citys, 1);
D = zeros(n, n);
for i = 1:n
    for j = 1:n
        if i~=j       %不等于
            D(i,j) = sqrt(sum((citys(i,:) - citys(j,:)).^2));
        else
            D(i,j) = 1e-4;
        end
    end
end

%% 初始化参数
m = 40;   %蚂蚁数量
alpha = 1;   %信息素重要程度因子
beta = 5;    %启发函数重要程度因子
rho = 0.1;   %信息素挥发因子
Q = 1;        %常系数
Eta = 1./D;    %启发函数
Tau = ones(n,n); %信息素矩阵
Table = zeros(m,n); %路径记录表
iter_max = 200;    %最大迭代次数
Route_best = zeros(iter_max, n);  %各代最佳路径
Length_best = zeros(iter_max,1);   %各代最佳路径长度
Length_ave = zeros(iter_max,1);    %各代路径的平均长度
m_c = 50;
trace_new = zeros(iter_max,1);
%% 蒙特卡罗实验
for k = 1:m_c
    trace = [];
    iter = 1;    %迭代次数初值
%% 迭代寻找最佳路径
while iter <= iter_max
    %随机产生各个蚂蚁的起点城市
    start = zeros(m, 1);
    for i = 1:m
        temp = randperm(n);
        start(i) = temp(1);
    end
    Table(:,1) = start;
    %构建解空间
    citys_index = 1:n;
    %逐个蚂蚁路径选择
    for i = 1:m
        %逐个城市路径选择
        for j = 2:n
            tabu = Table(i,1:(j-1));   %已访问的城市集合（禁忌表）
            allow_index = ~ismember(citys_index,tabu);
            allow = citys_index(allow_index);  %待访问的城市集合
            P = allow;
            %计算城市间的转移概率
            for k = 1:length(allow)
                P(k) = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
            end
            P = P/sum(P);
            %轮盘赌法选择下一个访问城市
            Pc = cumsum(P);   %累计和
            target_index = find(Pc>=rand);
            target = allow(target_index(1));
            Table(i,j) = target;
        end
    end
    %计算各个蚂蚁的路径距离
    Length = zeros(m,1);
    for i = 1:m
        Route = Table(i,:);
        for j = 1:(n-1)
            Length(i) = Length(i) + D(Route(j),Route(j+1));
        end
        Length(i) = Length(i) + D(Route(n),Route(1));
    end
    %计算最短路径距离及平均距离
    if iter == 1
        [min_Length,min_index] = min(Length);
        Length_best(iter) = min_Length;
        Length_ave(iter) = mean(Length);
        Route_best(iter,:) = Table(min_index,:);
    else
        [min_Length,min_index] = min(Length);
        Length_best(iter) = min(Length_best(iter - 1),min_Length);
        Length_ave(iter) = mean(Length);
        if Length_best(iter) == min_Length
            Route_best(iter,:) = Table(min_index,:);
        else
            Route_best(iter,:) = Route_best((iter-1),:);
        end
    end
    %更新信息素
    Delta_Tau = zeros(n,n);
    %逐个蚂蚁计算
    for i = 1:m
        %逐个城市计算
        for j = 1:(n-1)
            Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
        end
        Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
    end
    Tau = (1- rho)*Tau + Delta_Tau;
    %迭代次数加1，清空路径记录表
    trace = [trace;Length_best(iter)];
    iter = iter + 1;
    Table = zeros(m,n);
end
trace_new = trace_new + trace;
end
toc
%% 优化过程
figure;
plot(1:iter_max,trace_new./m_c);
title('蚁群算法50次蒙特卡罗收敛曲线')
xlabel('代数')
ylabel('最优值/km')
axis([0,200,14000,40000]);
grid on

%% 结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);

%% 绘图
% figure(1)
% plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
%      [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'o-');
%  grid on;
%  for i = 1:size(citys,1)
%      text(citys(i,1),citys(i,2),['  ' num2str(i)]);
%  end
% text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'    起点');
% text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'    终点');
% xlabel('城市位置横坐标');
% ylabel('城市位置纵坐标');
% title([ '蚁群算法优化路径（最短距离）：' num2str(Shortest_Length) ')' ]);
% % figure(2);
% % plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r');
% % legend('最短路径','平均距离');
% % xlabel('迭代次数');
% % ylabel('距离');
% % title('各代最短距离与平均距离对比');

%% 画出最优解的路线图
n = size(citys,1);
figure
hold on
plot([citys(Shortest_Route(1),1),citys(Shortest_Route(n),1)],[citys(Shortest_Route(1),2),...
    citys(Shortest_Route(n),2)],'ks-','Markersize',8,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
text(citys(1,1),citys(1,2),['   ' num2str(1)]);
hold on
for i=2:n
    plot([citys(Shortest_Route(i-1),1),citys(Shortest_Route(i),1)],[citys(Shortest_Route(i-1),2),...
        citys(Shortest_Route(i),2)],'ks-','Markersize',8,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
    text(citys(i,1),citys(i,2),['  ' num2str(i)]);
    hold on
end
legend('规划路径')
title([ '蚁群算法优化路径（最短距离）：' num2str(Shortest_Length)],'fontsize',10);
xlabel('城市位置横坐标/km','fontsize',10)
ylabel('城市位置纵坐标/km','fontsize',10)
grid on
x = citys(:,1);
y = citys(:,2);
startx=x(Shortest_Route(1)); %起点x坐标
starty=y(Shortest_Route(1)); %起点y坐标
endx=x(Shortest_Route(n));
endy=y(Shortest_Route(n));
text(startx,starty,'    起点'); %标记起点
text(endx,endy,'    终点')%标记终点
set(gca,'LineWidth',1.5);  %边框加粗,美观
axis([1000 1.1*max(x) 500 1.1*max(y)]); %设置尺寸大小




















