function [gbestval, gbestlog, n1, n2, state] = MHS_QPSO(ifu,n,c,MaxEFs,bu,bd)
BU = repmat(bu,1,c);
BD = repmat(bd,1,c);
%% Build the initial database
POP = lhsdesign(n,c).*(ones(n,1)*(BU-BD))+ones(n,1)*BD;
obj = compute_objectives(POP,c,ifu);
Data = [POP,obj];
gbestval = min(Data(:,c+1));
gbestlog = zeros(1,MaxEFs); 
startEFs = n;
currentEFs = startEFs;
gbestlog(1,1:startEFs) = repmat(gbestval,1,startEFs);
if c < 50
    gmax = max(c,20);
    pn = c;
else
    gmax = 50;
    pn = 50;
end
state = zeros(1,MaxEFs); % Search method record
a = [];
n1 = 0;n2 = 0;
%% Main Loop
while currentEFs < MaxEFs
    %% Search behavior selection scheme
    P = (currentEFs-startEFs)/(MaxEFs-startEFs);
    r = rand();
    if r > P
        %% Multi-region global search
        [ cbest ] = GS(Data,pn,c,bu,bd,gmax);
        %% Screen candidate solutions
        x = Data(:,1:c);
        y = Data(:,c+1);
        [ct] = reliability(x,y,cbest(:,1:c)); % Get the reliability of candidate solutions
        [~, index1] = min(cbest(:,c+1)); % Select the solution with the lowest predicted value
        [~, index2] = max(ct); % Select the solution with the highest reliability
        %% Real function evaluation
        cbest(index1,c+1) = compute_objectives(cbest(index1,1:c),c,ifu);
        currentEFs = currentEFs+1;
        Data = [Data;cbest(index1,:)];
        gbestlog(currentEFs) = min(Data(:,c+1));
        state(currentEFs) = 1;
        if index1 ~= index2
            cbest(index2,c+1) = compute_objectives(cbest(index2,1:c),c,ifu);
            currentEFs = currentEFs+1;
            Data = [Data;cbest(index2,:)];
            gbestlog(currentEFs) = min(Data(:,c+1));
        end
        n1 = n1+1;
    else
        %% Local search with dynamic boundary adjustment
        [ lbest ] = LS(Data,pn,c,bu,bd,gmax);
        %% Real function evaluation
        lbest(:,c+1) = compute_objectives(lbest(:,1:c),c,ifu);
        currentEFs = currentEFs+1;
        Data = [Data;lbest];
        gbestlog(currentEFs) = min(Data(:,c+1));
        state(currentEFs) = 2;
        n2 = n2+1;
    end
    gbestval = min(Data(:,c+1));
end
gbestlog=gbestlog(1:MaxEFs);
state=state(1:MaxEFs);
end

