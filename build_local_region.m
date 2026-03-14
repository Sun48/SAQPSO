function [localBD, localBU] = build_local_region(x, bu, bd)
% 修复了索引报错并融合了所有逻辑的版本

[sampleCount, dim] = size(x);

% --- 健壮性检查：确保 bu 和 bd 是 1 x dim 的行向量 ---
if numel(bu) == 1, bu = repmat(bu, 1, dim); end
if numel(bd) == 1, bd = repmat(bd, 1, dim); end
if size(bu, 1) > 1, bu = bu'; end % 转为行向量
if size(bd, 1) > 1, bd = bd'; end

globalSpan = bu - bd;

% 初始化：默认为精英样本的包络范围，但不超出全局边界
localBD = max(min(x, [], 1), bd); 
localBU = min(max(x, [], 1), bu);

if sampleCount == 0
    return;
end

%% A. 获取区间并排序 (你的原始片段逻辑)
eliteCount = min(sampleCount, max(10, min(50, sampleCount)));
elite = x(1:eliteCount, :);
lu = [min(elite, [], 1); max(elite, [], 1)]'; % D行2列

[sorted_intervals, sorted_index] = sortrows(lu, 1);
sorted_intervals(:, 3) = sorted_index; % 记录维度原始索引

%% B. 寻找公共交集区间 (你的核心算法)
dims_index = zeros(dim, dim);
common_intervals = zeros(0, 2);
row_index = 1;
column_index = 1;
current_interval = sorted_intervals(1, 1:2);
dims_index(row_index, column_index) = sorted_intervals(1, 3);

for i = 2:size(sorted_intervals, 1)
    next_interval = sorted_intervals(i, 1:2);
    if next_interval(1) <= current_interval(2)
        % 更新交集
        current_interval(1) = max(current_interval(1), next_interval(1));
        current_interval(2) = min(current_interval(2), next_interval(2));
        column_index = column_index + 1;
        dims_index(row_index, column_index) = sorted_intervals(i, 3);
    else
        % 保存旧区间，开启新区间
        if current_interval(2) >= current_interval(1)
            common_intervals = [common_intervals; current_interval];
        end
        row_index = row_index + 1;
        column_index = 1;
        dims_index(row_index, column_index) = sorted_intervals(i, 3);
        current_interval = next_interval;
    end
end
% 处理最后一组
if current_interval(2) >= current_interval(1)
    common_intervals = [common_intervals; current_interval];
end

%% C. 映射与保护带计算 (增强版保护逻辑)
valid_rows = find(dims_index(:, 1) ~= 0);
for i = 1:numel(valid_rows)
    % 找到这一组交集对应的所有维度索引
    value_set = dims_index(valid_rows(i), dims_index(valid_rows(i), :) ~= 0);
    interval = common_intervals(i, :);
    
    % 计算 Gap，防止空间塌陷
    intervalSpan = interval(2) - interval(1);
    % 针对这组维度的全局跨度
    currentGlobalSpan = globalSpan(value_set);
    
    gap = max(0.05 * currentGlobalSpan, 0.1 * intervalSpan);
    
    % 应用新的边界并确保不超全局界
    lower = max(bd(value_set), interval(1) - gap);
    upper = min(bu(value_set), interval(2) + gap);
    
    % 映射回 localBD/localBU
    localBD(value_set) = lower;
    localBU(value_set) = upper;
end

%% D. 最终兜底：防止任何维度出现 upper <= lower
degenerate = localBU <= localBD;
if any(degenerate)
    localBD(degenerate) = bd(degenerate);
    localBU(degenerate) = bu(degenerate);
end

end