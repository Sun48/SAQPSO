function [ gbest_log ] = GS(Data,n,c,bu,bd,gmax)
if c < 50
    len = 5*c;
else
    len = 100;
end
%% Obtain center points of regions
Nc = 5; % the number of clusters
gbest_log = []; 
[~, center] = kmeans(Data(:,1:c),Nc); 
X = Data(:,1:c);
Y = Data(:,c+1);
dis = zeros(1,Nc-1); 
%% Multi-region search
for i = 1:Nc
    %% Obtain samples of model building
    ct = center;
    ct(i,:) = [];
    for j = 1:Nc-1
        dis(j) = norm((center(i,:)-ct(j,:)));
    end
    r = 0.5*max(dis)/((c^0.5)*(Nc-1)^(1/c));
    BD = (center(i,:)-r);
    BU = (center(i,:)+r);
    for j = 1:c
        if BD(j) < bd
            BD(j) = bd;
        end
        if BU(j) > bu
            BU(j) = bu;
        end
    end
    d = pdist2(center(i,:),X)';
    [~, idx] = sort(d);
    X = X(idx,:);
    Y = Y(idx,:);
    tt = 0;
    for j = 1:c
        tt = tt+(X(len,j)<BU(j)&&X(len,j)>BD(j));
    end
    if tt < c
        x = X(1:len,:);
        y = Y(1:len,:);
    else
        for k = len:size(X,1)
            tt = 0;
            for j = 1:c
                tt = tt+(X(k,j)<BU(j)&&X(k,j)>BD(j));
            end
            if tt < c
                break;
            end
        end
        x = X(1:k,:);
        y = Y(1:k,:);
    end 
    %% RBF
    srgtOPTRBF  = srgtsRBFSetOptions(x, y);
    srgtSRGTRBF = srgtsRBFFit(srgtOPTRBF);
    %% Obtain the initial particles
    pp = lhsdesign(n,c).*(ones(n,1)*(BU-BD))+ones(n,1)*BD; % Initialize particles
    fitness = srgtsRBFEvaluate(pp, srgtSRGTRBF);
    p = [pp, fitness];
    lbest = p;
    [~, Ib] = min(p(:,c+1));
    gbest = p(Ib,:);
    g = 0;
    %% QPSO
    while g < gmax
        [ p ] = move_qpso(p,BU,BD,gbest,lbest,g,gmax);
        fitness = srgtsRBFEvaluate(p(:,1:c), srgtSRGTRBF);
        p(:,c+1) = fitness;
        [best, Ib] = min(p(:,end));
        if best <= gbest(end)
            gbest = p(Ib,:);
        end
        I = find(p(:,end)<=lbest(:,end));
        lbest(I,:) = p(I,:);
        g = g+1;
    end
    gbest_log = [gbest_log; gbest];
end
end

