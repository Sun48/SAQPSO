function [ gbest ] = LS(Data, n, c, bu, bd, gmax)
% Data: 输入样本矩阵，最后一列为适应度
% n: 粒子数量
% c: 维度数量
% bu/bd: 全局上下界
% gmax: 最大迭代次数

xmax_t = zeros(1,c); 
xmin_t = zeros(1,c); 

%% 1. 获取模型构建样本
P = sortrows(Data, c+1); % 按适应度排序
% 计算到当前最优点的距离，筛选局部样本
dis = sqrt(sum((P(:,1:c)-P(1,1:c)).^2, 2));
[~, i] = sort(dis);
np = max(2, floor(0.2*size(Data,1)));
Datal = P(i(1:np),:);
x = Datal(:,1:c);
y = Datal(:,c+1);

%% 2. 拟合 RBF 代理模型
srgtOPTRBF  = srgtsRBFSetOptions(x, y);
srgtSRGTRBF = srgtsRBFFit(srgtOPTRBF);

%% 3. 构建局部有前途区域 (调用整合后的函数)
[localBD, localBU] = build_local_region(x, bu, bd);

%% 4. 初始化 QPSO 粒子
nn = min(np, n);
p = Datal(1:nn, :); % 取精英点作为初始粒子
lbest = p;
[~, Ib] = min(p(:, c+1));
gbest = p(Ib, :);
g = 0;

%% 5. QPSO 优化循环
while g < gmax
    % 动态边界调整：结合局部区域和粒子分布
    for t = 1:c
        xmax_t(t) = max(p(:,t));
        xmin_t(t) = min(p(:,t));
        span_t = xmax_t(t)-xmin_t(t);
        
        % 适度外扩粒子分布范围
        xmax_temp = xmax_t(t) + 0.25*span_t;
        xmin_temp = xmin_t(t) - 0.25*span_t;
        
        % 严格限制在 build_local_region 计算的范围内
        xmax_t(t) = min(localBU(t), xmax_temp);
        xmin_t(t) = max(localBD(t), xmin_temp);
        
        % 退化保护
        if xmax_t(t) <= xmin_t(t)
            xmax_t(t) = localBU(t);
            xmin_t(t) = localBD(t);
        end
    end
    
    BU = xmax_t;
    BD = xmin_t;
    
    % 粒子移动与评估
    [ p ] = move_qpso(p, BU, BD, gbest, lbest, g, gmax);
    f = srgtsRBFEvaluate(p(:,1:c), srgtSRGTRBF);
    p(:, c+1) = f;
    
    % 更新个体最好和全局最好
    [best, Ib] = min(p(:, end));
    if best <= gbest(end)
        gbest = p(Ib, :);
    end
    I = find(p(:, end) <= lbest(:, end));
    lbest(I, :) = p(I, :);
    
    g = g + 1;
end
end


% lu = [min(lhx); max(lhx)]; % 获取局部空间的上下界（计算局部有前途区域）,对这50个X以列为单位，确定每一列的上界和下界，生成的是2行D列的矩阵
% [sorted_intervals,sorted_index] = sortrows(lu', 1); % 按照升序对局部空间的边界进行排序，按照下界从小到大排
% common_intervals = []; % 初始化公共区间，用于存储最终的核心空间区间
% dims_index=zeros(Dim,Dim); % 初始化维度索引矩阵 % 标志位
% current_interval = sorted_intervals(1, 1:2); % 正在处理的空间，获取排序后的sorted_intervals第一行的两个数值（最小的下界+对应的上界）
% sorted_intervals(:,3)=sorted_index; % 给排序后的区间添加索引,添加到sorted_intervals的第三列
% %% 然后，判断当前空间是否和下一个空间有交集，直到当前空间和下一个空间没有交集（当前空间不断较小）其中，第一次的当前空间取排序后的第一个即空间下界最小的
% %% 表面是当前空间只跟下一个空间交，其实他们会将交集作为当前空间，整合后再去跟下一个交，也相当于所有的都交了一下
% row_index=1; % 初始化行索引，表示当前维度所在的行
% conclumn_index=1; % 初始化列索引，表示当前维度所在的列
% dims_index(row_index,conclumn_index)=sorted_intervals(1,3); % 将第一个排序后的区间索引存储到 dims_index 矩阵的第一个位置
% for i = 2:size(sorted_intervals, 1) % 循环处理每一个排序后的区间，从第二个开始
% next_interval = sorted_intervals(i, 1:2); % 获取下一个区间的上下界（即下一个维度的上下界）
% if next_interval(1) <= current_interval(2) % 如果当前区间和下一个区间有交集，更新交集区间
% current_interval(1) = max(current_interval(1), next_interval(1)); % 更新当前区间的下界为两个区间下界的较大值
% current_interval(2) = min(current_interval(2), next_interval(2)); % 更新当前区间的上界为两个区间上界的较小值
% conclumn_index=conclumn_index+1; % 列索引递增，表示当前维度索引列向右移动
% dims_index(row_index,conclumn_index)=sorted_intervals(i,3); % 将当前维度的索引存入 dims_index 矩阵
% else % 如果当前区间和下一个区间没有交集，表示需要划分新的一行
% conclumn_index=1;
% row_index=row_index+1;
% dims_index(row_index,conclumn_index)=sorted_intervals(i,3); % 将新维度的索引存入新的行列中
% if current_interval(2) >= current_interval(1) % 如果当前区间有效，将其加入到公共区间中,就是说这个区间没啥问题，上界要大于等于下界
% common_intervals = [common_intervals; current_interval]; % 将当前的交集区间存储到公共区间数组中
% end
% current_interval = next_interval; % 更新当前区间为下一个区间
% %disp(['娌℃湁浜ら泦锛?=',num2str(i)]);
% end
% end
% %% 最后，将所有符合条件的维度存储到公共索引集中，然后为了防止空间过小加了个gap
% if current_interval(2) >= current_interval(1)
% common_intervals = [common_intervals; current_interval];
% end
% K=0; % 初始化
% raws=length(find(dims_index(:,1)~=0)); % 获取维度索引的非零行数,返回第一列非零元素的个数，我感觉是控制搜索行数
% for i=1:raws
% value_set=dims_index(i,find(dims_index(i,:)~=0)); % 获取该行的非零值
% conclumns=length(value_set); % 计算该行的列数，控制后面几个列是一样的上下界
% if conclumns<K