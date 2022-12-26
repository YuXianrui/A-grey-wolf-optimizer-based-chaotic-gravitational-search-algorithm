function G = Gconstant(iteration,max_it)
%GC 此处显示有关此函数的摘要
%   此处显示详细说明
alfa = 20; G0=100;
G=G0*exp(-alfa*iteration/max_it);
end

