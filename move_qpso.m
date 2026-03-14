function [ POP ] = move_qpso(POP,VRmax,VRmin,g_best,lbest,i,me)
ps=size(POP,1);
D=size(VRmax,2);
pos=POP(:,1:D); 

pbest=lbest(:,1:D);
gbest=g_best(:,1:D);
pavg=mean(pbest); 
pavgrep=repmat(pavg,ps,1);
gbestrep=repmat(gbest,ps,1);

fi=rand(ps,D);
p=fi.*pbest+(1-fi).*gbestrep;
alpha=(1-0.5)*(me-i)/me+0.5;

b=alpha*abs(pavgrep-pos);
u=rand(ps,D);
pos=p+((-1).^ceil(0.5+rand(ps,D))).*b.*(-log(u)); 
pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
    +((pos<VRmin)|(pos>VRmax)).*(VRmin+(VRmax-VRmin).*rand(ps,D)); 
POP(:,1:D)=pos;
end