function  [ y ]=reliability(Xtrain,Ytrain,Xtest)
[m,~]=size(Xtest);
K=4; % number of adjacent evaluation points
for i=1:m
    d=pdist2(Xtest(i,:),Xtrain);
    [val,ind]=sort(d);
    dmin(i)=val(1);
    dg(i)=mean(val(1:K));
    t1=Ytrain(ind(1:K));
    va(i)=sqrt((var(t1)));
    
    a=Xtest(i,:);
    norm_p1 = norm(a);
    for j=1:K
        b=Xtrain(ind(j),:);
        dot_product = dot(a,b);
        norm_p2 = norm(b);
        cos_theta(j) = dot_product / (norm_p1 * norm_p2);
    end
    cos(i)=mean(cos_theta);
end
dmin=dmin/sum(dmin);
dg=dg/sum(dg);
va=va/sum(va);
cos=cos/sum(cos);
y=-dmin.*(va+dg)./cos;
end

