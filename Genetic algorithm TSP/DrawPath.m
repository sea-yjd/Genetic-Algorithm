function DrawPath(Chrom,X,min)
%% 画路线图函数
% 输入
% Chrom     待画路线
% X    各城市的坐标位置

R = [Chrom(1,:) Chrom(1,1)];   %一个随机解（个体）
figure;
hold on
plot(X(:,1), X(:,2),'o','color',[0.5,0.5,0.5])
plot(X(Chrom(1,1),1),X(Chrom(1,1),2),'rv','MarkerSize',20)
for i = 1:size(X,1)
    text(X(i,1)+0.05,X(i,2)+0.05,num2str(i),'color',[1,0,0]);
end
A = X(R,:);
row = size(A,1);
for i = 2:row
    [arrowx,arrowy] = dsxy2figxy(gca,A(i-1:i,1),A(i-1:i,2));   %坐标转换
    annotation('textarrow',arrowx,arrowy,'HeadWidth',8,'color',[0,0,1]);
end
hold off
xlabel('城市位置横坐标');
ylabel('城市位置纵坐标');
title([ '遗传算法优化路径（最短距离）：' num2str(min)]);
box on





