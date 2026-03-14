function [ bu, bd ] = func_boundary(ifu)
if ifu==1
    bu=5.12;bd=-5.12;
end
if ifu==2
    bu=2.048;bd=-2.048;
end
if ifu==3
    bu=32.768;bd=-32.768;
end
if ifu==4
    bu=600;bd=-600;
end
if ifu==5
    bu=5;bd=-5;
end
if ifu==6
    bu=5;bd=-5;
end
if ifu==7
    bu=5;bd=-5;
end
end

