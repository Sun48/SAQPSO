warning('off')
format compact;
addpath('CEC2005');
clear all;
clc;
func_num=7; 
runs=5; 
%% Record
ftime=zeros(func_num,runs);
ft_mean=zeros(1,func_num); 
fbest=zeros(func_num,runs);
fb_mean=zeros(1,func_num); 
fb_best=zeros(1,func_num); 
fb_std=zeros(1,func_num); 
for c = [30] % Dimensions
        if c < 50
            n=5*c;
            MaxEFs=11*c;
        else
            n=100;
            MaxEFs=1000;
        end
        fb_log=zeros(runs,MaxEFs);
        state_log=zeros(runs,MaxEFs);
        fb_log_mean=zeros(func_num,MaxEFs);
        %% Main Loop
        for i = 1:func_num
            [bu,bd]=func_boundary(i);
            for j = 1:runs
                T=cputime;
                [gbestval,gbestlog,n1,n2,state]=MHS_QPSO(i,n,c,MaxEFs,bu,bd);
                ftime(i,j)=cputime-T; 
                fbest(i,j)=gbestval;
                fb_log(j,:)=gbestlog;
                state_log(j,:)=state;
                fprintf('**********func_num=%d, runs=%d, best_value=%g, GS=%g, LS=%g**********\n\n',i,j,gbestval,n1,n2);
            end
            ft_mean(i)=mean(ftime(i,:));
            fb_mean(i)=mean(fbest(i,:));
            fb_best(i)=min(fbest(i,:));
            fb_std(i)=std(fbest(i,:)); 
            fb_log_mean(i,:)=mean(fb_log,1);
        end
        %% Save Result（修正版：7行5列，对齐每个函数的统计结果）
        % 1. 确保result_test目录存在
        if ~exist('result_test', 'dir')
            mkdir('result_test');
        end

        % 2. 生成Excel文件名（兼容所有MATLAB版本）
        time_str = datestr(now, 'yyyy-mm-dd-HH-MM');
        excel_file = fullfile('result_test', ['result_', num2str(c), 'd_', time_str, '.xlsx']);

        % 3. 核心数据：7行5列（每行对应1个函数，每列对应1个指标）
        core_data = [
            (1:func_num)', ...   % 第1列：函数编号（7×1）
            ft_mean', ...        % 第2列：平均耗时（7×1）
            fb_mean', ...        % 第3列：平均最优值（7×1）
            fb_best', ...        % 第4列：全局最优值（7×1）
            fb_std' ...          % 第5列：标准差（7×1）
            ];  % 注意：用逗号`,`横向拼接，最终是7行5列的矩阵

        % 表头
        core_header = {'Func_Num', 'FT_Mean', 'FB_Mean', 'FB_Best', 'FB_Std'};

        % 4. 写入Excel：表头在A1，数据从A2开始（刚好7行，对应7个函数）
        xlswrite(excel_file, core_header, '核心结果', 'A1');  % 写表头
        xlswrite(excel_file, core_data, '核心结果', 'A2');    % 写7行数据

        % （可选）保存其他详细数据到不同Sheet
        xlswrite(excel_file, ftime, '各Run耗时');
        xlswrite(excel_file, fbest, '各Run最优值');
        xlswrite(excel_file, fb_log_mean, '迭代日志均值');

        fprintf('维度 %dd 结果已保存至：%s\n', c, excel_file);
end


