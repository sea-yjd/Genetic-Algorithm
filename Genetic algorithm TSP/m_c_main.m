clear
clc
close all
load citys_data.mat
tic

NIND = 100;   %种群大小
MAXGEN = 200;
Pc = 0.75;    %交叉概率
Pm = 0.01;   %变异概率
GGAP = 0.8;    %代沟
D = Distanse(citys);   %生成距离矩阵
N = size(D,1);     %(34*34)
m_c = 50;
trace_new = zeros(MAXGEN,1);
%% 蒙特卡罗实验
for m = 1:m_c
%% 初始化种群
Chrom = InitPop(NIND,N);
trace = [];
%% 在二维图上画出所有坐标点
% figure 
% plot(citys(:,1),citys(:,2),'o');
%% 画出随机解的路线图
% DrawPath(Chrom(1,:),citys);
% pause(0.0001);
%% 输出随机解的路线和总距离
% disp('初始种群中的一个随机值：');
% OutputPath(Chrom(1,:));
% Rlength = PathLength(D, Chrom(1,:));
% disp(['总距离：',num2str(Rlength)]);
% disp('================================');
%% 优化
gen = 0;
% figure;
% hold on;box on
% xlim([0,MAXGEN]);
% title('优化过程');
% xlabel('代数');
% ylabel('最优值');
ObjV = PathLength(D,Chrom);   %计算路线长度
preObjV = min(ObjV);
while gen < MAXGEN
    %% 计算适应度
    ObjV = PathLength(D, Chrom);   %计算路线长度
    %fprint('%d    %1.10f\n',gen, min(ObjV))
    trace = [trace;min(ObjV)];
%     line([gen-1,gen],[preObjV,min(ObjV)]);
%     pause(0.0001)
    preObjV = min(ObjV);
    FitnV = Fitness(ObjV);
    %% 选择
    SelCh = Select(Chrom, FitnV,GGAP);
    %% 交叉操作
    SelCh = Recombin(SelCh,Pc);
    %% 变异
    SelCh = Mutate(SelCh,Pm);
    %% 逆转操作
    SelCh = Reverse(SelCh, D);
    %% 重插入子代的新种群
    Chrom = Reins(Chrom,SelCh,ObjV);
    
    %% 更新迭代次数
    gen = gen+1;
end
trace_new = trace_new + trace;
end
toc


%% 画出最优解的优化过程
figure
plot(1:MAXGEN,trace_new./m_c);      %画出每代的最优点
xlabel('进化代数');
ylabel('最优值/km');
title('遗传算法50次蒙特卡罗收敛曲线');
grid on;
%% 画出最优解的路线图
ObjV = PathLength(D,Chrom);    %计算路线长度
[minObjV,minInd] = min(ObjV);
% DrawPath(Chrom(minInd(1),:),citys,ObjV(minInd(1)));
tourGbest = Chrom(minInd,:);
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
title([ '遗传算法优化路径（最短距离）：' num2str(ObjV(minInd(1)))],'fontsize',10);
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
p = OutputPath(Chrom(minInd(1),:));
disp(['总距离：',num2str(ObjV(minInd(1)))]);
disp('==============================');















 
