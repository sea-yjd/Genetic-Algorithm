function D = Distanse(a)
%% 计算两两城市间距离
% 输入 a  各城市的位置坐标
% 输出D    两两城市间距离

row = size(a,1);
D = zeros(row,row);
for i = 1:row
    for j = i+1:row
        D(i,j) = ((a(i,1)-a(j,1))^2 + (a(i,2)-a(j,2))^2)^0.5;
        D(j,i) = D(i,j);
    end
end