function [ obj ] =compute_objectives(POP,c,problem_num)
obj=[];
global initial_flag
initial_flag=0;
switch problem_num
    case 1 %'Ellipsoid' Bounds[-5.12,5.12]  f_bias=0
        P=[1:c];P=ones(size(POP,1),1)*P;
        obj=sum(P.*(POP(:,1:c).^2),2);     
    case 2 %'Rosenbrock' Bounds[-2.048,2.048]  f_bias=0
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        obj=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2);
    case 3 %'Ackley' Bounds[-32.768,32.768]  f_bias=0
        a=20;b=0.2;cc=2*pi;
        obj=0-a*exp(0-b*(mean(POP(:,1:c).^2,2)).^0.5)-exp(mean(cos(cc*POP(:,1:c)),2))+a+exp(1);
    case 4 %'Griewank' Bounds[-600,600]  f_bias=0
        P=[1:c].^0.5;P=ones(size(POP,1),1)*P;
        obj=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
    case 5 %'Rastrigin' Bounds[-5,5]  f_bias=0
        obj=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
    case 6 % 10.Shifted Rotated Rastrigin's  Function  Bounds[-5,5]  f_bias=-330
        obj=benchmark_func(POP(:,1:c),10);
    case 7 % 19.Rotated Hybrid Composition Function 2 with a Narrow Basin for the Global Optimum  Bounds[-5,5]]  f_bias=10 
        obj=benchmark_func(POP(:,1:c),19);
end
end