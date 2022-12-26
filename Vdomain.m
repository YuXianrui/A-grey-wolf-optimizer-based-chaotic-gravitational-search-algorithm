%Location update strategy and preconditions for execution
function X = Vdomain(M,X,BestChart,iteration,MeanChart,max_it, ac ,V)
[N,dim]=size(X);
[Ms,ds]=sort(M,'descend');
pu=3; %Steepness coefficient
u = 1/exp(1);
a=2-iteration*((2)/max_it); % a decreases linearly fron 2 to 0
Epsilon = 0.1;
if iteration >=0.5*max_it % n is the number of affected agents
    n = ceil(N * ((iteration/(pu*max_it)+ u)*log(((iteration/(pu*max_it)+ u)) - u * log(u))));
else
    n = floor(N * ((iteration/(pu*max_it)+ u)*log(((iteration/(pu*max_it)+ u)) - u * log(u))));
end

if iteration >=2
    kvalue1=BestChart(1,iteration);
    kvalue2=BestChart(1,iteration-1);
    kvalue_1=kvalue2-kvalue1;
    OP =kvalue_1/(MeanChart(1,iteration)-kvalue1);
    if OP < Epsilon %Preconditions for execution
        for i = 1:N %Location update strategy addressed in GWO
            if i>=1 && i<=n
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
                
                A1=2*a*r1-a; % Equation (3.3)
                C1=2*r2; % Equation (3.4)
                
                D_alpha=abs(C1*(X(ds(1),:))-X(ds(N-i+1),:)); % Equation (3.5)-part 1
                X1=(X(ds(1),:))-A1*D_alpha; % Equation (3.6)-part 1
                
                r1=rand();
                r2=rand();
                
                A2=2*a*r1-a; % Equation (3.3)
                C2=2*r2; % Equation (3.4)
                
                D_beta=abs(C2*(X(ds(2),:))-X(ds(N-i+1),:)); % Equation (3.5)-part 2
                X2=(X(ds(2),:))-A2*D_beta; % Equation (3.6)-part 2
                
                r1=rand();
                r2=rand();
                
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)
                
                D_delta=abs(C3*(X(ds(3),:))-X(ds(N-i+1),:)); % Equation (3.5)-part 3
                X3=(X(ds(3),:))-A3*D_delta; % Equation (3.5)-part 3
                
                X(ds(N-i+1),:)=(X1+X2+X3)/3;% Equation (3.7)
            else
                V=rand(N,dim).*V+ac;
                X(ds(N-i+1),:) =X(ds(N-i+1),:)+V(ds(N-i+1),:);
            end
           
        end
        
    else
        
        V=rand(N,dim).*V+ac;
        X=X+V;
        
    end
end


