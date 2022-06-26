clear
clc
close all
load citys_data.mat
tic
cityCoor = [citys(:,1) citys(:,2)];     %城市坐标
% figure
% plot(cityCoor(:,1),cityCoor(:,2),'ms','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g');
% legend('城市位置');
% title('城市分布图','fontsize',12)
% xlabel('km','fontsize',12)
% ylabel('km','fontsize',12)
% ylim([min(cityCoor(:,2))-1 max(cityCoor(:,2))+1])
% grid on


%% 计算城市间距离
n = size(cityCoor,1);    %城市数目
cityDist = zeros(n,n);    %城市距离矩阵
for i = 1:n
    for j = 1:n
        if i~=j
            cityDist(i,j) = ((cityCoor(i,1)-cityCoor(j,1))^2 + ...
                (cityCoor(i,2)-cityCoor(j,2))^2)^0.5;
        end
        cityDist(j,i) = cityDist(i,j);
    end
end

nMax = 200;   %进化次数
indiNumber = 100;     %个体数目
individual = zeros(indiNumber,n);

%% 粒子初始化
for i = 1:indiNumber
    individual(i,:) = randperm(n);    %粒子位置
end

%% 计算种群适应度
indiFit = fitness(individual,cityCoor, cityDist);
[value, index] = min(indiFit);
tourPbest = individual;      %当前个体最优
tourGbest = individual(index,:);     %当前全局最优
recordPbest = inf * ones(1,indiNumber);    %个体最优纪录
recordGbest = indiFit(index);    %群体最优纪录
xnew1 = individual;

%% 循环寻找最优路径
L_best=zeros(1,nMax);
for N=1:nMax    %N表示迭代次数
    %计算适应度值
    indiFit=fitness(individual,cityCoor,cityDist);
    
    %更新当前最优和历史最优
    for i=1:indiNumber
        if indiFit(i)<recordPbest(i)
            recordPbest(i)=indiFit(i);
            tourPbest(i,:)=individual(i,:);
        end
        if indiFit(i)<recordGbest
            recordGbest=indiFit(i);
            tourGbest=individual(i,:);
        end
    end
    
    [value,index]=min(recordPbest);
    recordGbest(N)=recordPbest(index);
    %% 交叉操作
    for i = 1:indiNumber
        %% 与个体最优进行交叉
        c1 = unidrnd(n-1);      %产生交叉位
        c2 = unidrnd(n-1);      %产生交叉位
        while c1 == c2
            c1 = round(rand * (n-2)) + 1;
            c2 = round(rand * (n-2)) + 1;
        end
        chb1 = min(c1,c2);
        chb2 = max(c1,c2);
        cros = tourPbest(i,chb1:chb2);    %交叉区域矩阵
        ncros = size(cros,2);       %交叉区域元素个数
        %删除与交叉区域相同的元素
        for j = 1:ncros
            for k = 1:n
                if xnew1(i,k) == cros(j)
                    xnew1(i,k) = 0;
                    for t = 1:n-k
                        temp = xnew1(i,k+t-1);
                        xnew1(i,k+t-1) = xnew1(i,k+t);
                        xnew1(i,k+t) = temp;
                    end
                end
            end
        end
        %插入交叉区域
        xnew1(i,n-ncros+1:n) = cros;
        %新路径长度短则接受
        dist = 0;
        for j = 1:n-1
            dist = dist + cityDist(xnew1(i,j),xnew1(i,j+1));
        end
        dist = dist + cityDist(xnew1(i,1),xnew1(i,n));
        if indiFit(i) > dist
            individual(i,:) = xnew1(i,:);
        end
        % 与全体最优进行交叉
        c1=round(rand*(n-2))+1;  %产生交叉位
        c2=round(rand*(n-2))+1;  %产生交叉位
        while c1==c2
            c1=round(rand*(n-2))+1;
            c2=round(rand*(n-2))+1;
        end
        chb1=min(c1,c2);
        chb2=max(c1,c2);
        cros=tourGbest(chb1:chb2); 
        ncros=size(cros,2);      
        %删除与交叉区域相同元素
        for j=1:ncros
            for k=1:n
                if xnew1(i,k)==cros(j)
                    xnew1(i,k)=0;
                    for t=1:n-k
                        temp=xnew1(i,k+t-1);
                        xnew1(i,k+t-1)=xnew1(i,k+t);
                        xnew1(i,k+t)=temp;
                    end
                end
            end
        end
        %插入交叉区域
        xnew1(i,n-ncros+1:n)=cros;
        %新路径长度变短则接受
        dist=0;
        for j=1:n-1
            dist=dist+cityDist(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+cityDist(xnew1(i,1),xnew1(i,n));
        if indiFit(i)>dist
            individual(i,:)=xnew1(i,:);
        end
    

        %% 变异操作
        c1 = round(rand * (n-1)) + 1;    %产生变异位
        c2 = round(rand * (n-1)) + 1;    %产生变异位
        while c1==c2
            c1 = round(rand * (n-2)) + 1;
            c2 = round(rand * (n-2)) + 1;
        end
        temp = xnew1(i,c1);
        xnew1(i,c1) = xnew1(i,c2);
        xnew1(i,c2) = temp;

        % 新路径长度变短则接受
        dist = 0;
        for j = 1:n-1
            dist = dist + cityDist(xnew1(i,j),xnew1(i,j+1));
        end
        dist = dist + cityDist(xnew1(i,1),xnew1(i,n));
        if indiFit(i) > dist
            individual(i,:) = xnew1(i,:);
        end
    end
    [value,index]=min(indiFit);
    L_best(N)=indiFit(index);
    tourGbest=individual(index,:); 
    
end
toc

%% 结果作图
figure
plot(L_best)
title('混合粒子群算法优化过程')
xlabel('代数')
ylabel('最优值/km')
axis([0,200,14000,40000]);
grid on


figure

hold on
plot([cityCoor(tourGbest(1),1),cityCoor(tourGbest(n),1)],[cityCoor(tourGbest(1),2),...
    cityCoor(tourGbest(n),2)],'ks-','Markersize',8,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
text(cityCoor(1,1),cityCoor(1,2),['   ' num2str(1)]);
hold on
for i=2:n
    plot([cityCoor(tourGbest(i-1),1),cityCoor(tourGbest(i),1)],[cityCoor(tourGbest(i-1),2),...
        cityCoor(tourGbest(i),2)],'ks-','Markersize',8,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
    text(cityCoor(i,1),cityCoor(i,2),['  ' num2str(i)]);
    hold on
end
legend('规划路径')
title([ '混合粒子群算法优化路径（最短距离）：' num2str(L_best(:,nMax))],'fontsize',10);
xlabel('城市位置横坐标/km','fontsize',10)
ylabel('城市位置纵坐标/km','fontsize',10)

grid on
x = citys(:,1);
y = citys(:,2);
disp(['最短距离:' num2str(L_best(:,nMax))]);
disp(['最短路径:' num2str( [tourGbest tourGbest(1)] )]);
startx=x(tourGbest(1)); %起点x坐标
starty=y(tourGbest(1)); %起点y坐标
endx=x(tourGbest(n));
endy=y(tourGbest(n));
text(startx,starty,'    起点'); %标记起点
text(endx,endy,'    终点')%标记终点
set(gca,'LineWidth',1.5);  %边框加粗,美观
axis([1000 1.1*max(x) 500 1.1*max(y)]); %设置尺寸大小

%% 画出最优解的路线图
% ObjV = PathLength(cityDist,tourGbest);    %计算路线长度
% [minObjV,minInd] = min(ObjV);
% DrawPath(tourGbest(minInd(1),:),citys,ObjV(minInd(1)));









