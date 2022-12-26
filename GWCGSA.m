function [Fbest, Lbest ,BestChart]=GWCGSA(N,max_it,F_index)
%V:   Velocity.
%a:   Acceleration.
%M:   Mass.  Ma=Mp=Mi=M;
%dim: Dimension of test function.
%N:   Number of agents.
%X:   Position of agents. dim-by-N matrix.
%R:   Distance between agents in search space.
%[low-up]: Allowable range for search space.
%Rnorm:  Norm in eq.8.
%Rpower: Power of R in eq.7.

ElitistCheck=1;
G_History=zeros(1,max_it);
Rnorm=2;
Rpower=1;
min_flag=1; % 1: minimization, 0: maximization
chValueInitial=20;
chaosIndex=10;
low=-100;
up=100;
dim=50;

%random initialization for agents.
X=initialization(dim,N,up,low);

%create chart of best so far and average fitnesses.
BestChart=[];MeanChart=[];

V=zeros(N,dim);
wMax=chValueInitial;
wMin=1e-10;
for iteration=1:max_it
    chValue=wMax-iteration*((wMax-wMin)/max_it);
    
    %     iteration
    
    %Checking allowable range.
    X=space_bound(X,up,low);
    
    %Evaluation of agents.
    fitness=cec14_func(X',F_index);
    
    [best best_X]=min(fitness); %min: minimization.  max: maximization.
    
    if iteration==1
        Fbest=best;Lbest=X(best_X,:);
    end
    if best<Fbest  % < : minimization. > : maximization
        Fbest=best;Lbest=X(best_X,:);
    end
    
    BestChart=[BestChart Fbest];
    MeanChart=[MeanChart mean(fitness)];
    
    %Calculation of M. eq.14-20
    [M]=massCalculation(fitness,min_flag);
    
    %Calculation of Gravitational constant. eq.13.
    G=Gconstant(iteration,max_it);

    switch chaosIndex
        case 1
            G=Gconstant(iteration,max_it);
        case 2
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 3
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 4
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 5
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 6
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 7
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 8
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 9
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 10
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
        case 11
            G=G+chaos(chaosIndex-1,iteration,max_it,chValue);
    end
    G_History(iteration)=G;
    test_G(iteration)=G;
    if iteration==499
        alisop=0;
    end
    %Calculation of accelaration in gravitational field. eq.7-10,21.
    [a,kbest]=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);
    
    %Agent movement
    [X,V]=move(X,a,V);
    
    %Location update strategy and preconditions for execution
    X = Vdomain(M,X,BestChart,iteration,MeanChart,max_it,a,V);
end %iteration



