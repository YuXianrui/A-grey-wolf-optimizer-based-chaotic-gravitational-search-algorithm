function G = Gconstant(iteration,max_it)
%GC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
alfa = 20; G0=100;
G=G0*exp(-alfa*iteration/max_it);
end

