%  Snow ablation optimizer (SAO)                                                                   
%                                                                                                     
%  Developed in MATLAB R2018a                                                                 
%                                                                                                     
%  Author : Lingyun Deng, Sanyang Liu                                                         
%                                                                                                     
%         e-Mail: lingyundeng@stu.xidian.edu.cn                                                            
%                 syliu@xidian.edu.cn 
                                                                                                                                                                                                        
%  Main paper:                                                                                        
%  Lingyun Deng,Sanyang Liu. Snow ablation optimizer: a novel metaheuristic technique for numerical optimization and engineering design
%  Journal name= Expert Systems with Applications, DOI: https://doi.org/10.1016/j.eswa.2023.120069
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of iterations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers

% To run SAO: [Best_pos,Best_score,Convergence_curve]=SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%______________________________________________________________________________________________


function [Best_pos,Best_score,Convergence_curve]=SAO(N,Max_iter,lb,ub,dim,fobj)

if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
%Initialize the set of random solutions
X=initialization_SAO(N,dim,ub,lb);

Best_pos=zeros(1,dim);
Best_score=inf;
Objective_values = zeros(1,size(X,1));

Convergence_curve=[];
N1=floor(N*0.5);
Elite_pool=[];

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    Objective_values(1,i)=fobj(X(i,:));
    if i==1
        Best_pos=X(i,:);
        Best_score=Objective_values(1,i);
    elseif Objective_values(1,i)<Best_score
        Best_pos=X(i,:);
        Best_score=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end

[~,idx1]=sort(Objective_values);
second_best=X(idx1(2),:);
third_best=X(idx1(3),:);
sum1=0;
for i=1:N1
    sum1=sum1+X(idx1(i),:);
end
half_best_mean=sum1/N1;
Elite_pool(1,:)=Best_pos;
Elite_pool(2,:)=second_best;
Elite_pool(3,:)=third_best;
Elite_pool(4,:)=half_best_mean;

Convergence_curve(1) = Best_score;

for i=1:N
    index(i)=i;
end

Na=N/2;
Nb=N/2;

%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while l<=Max_iter
    RB=randn(N,dim);          %Brownian random number vector
    T=exp(-l/Max_iter);
    k=1;
    DDF=0.35*(1+(5/7)*(exp(l/Max_iter)-1)^k/(exp(1)-1)^k);
    M=DDF*T;
    
    %% Calculate the centroid position of the entire population
    for j=1:dim
        sum1=0;
        for i=1:N
            sum1=sum1+X(i,j);
        end
        X_centroid(j)=sum1/N;
    end
    
    %% Select individuals randomly to construct pop1 and pop2
    index1=randperm(N,Na);
    index2=setdiff(index,index1);
    
    for i=1:Na
        r1=rand;
        k1=randperm(4,1);
        for j=1:size(X,2) % in j-th dimension
            X(index1(i),j)= Elite_pool(k1,j)+RB(index1(i),j)*(r1*(Best_pos(j)-X(index1(i),j))+(1-r1)*(X_centroid(j)-X(index1(i),j)));
        end
    end
    
    if Na<N
    Na=Na+1;
    Nb=Nb-1;
    end

    
    if Nb>=1
    for i=1:Nb
        r2=2*rand-1;
        for j=1:size(X,2) % in j-th dimension
            X(index2(i),j)= M*Best_pos(j)+RB(index2(i),j)*(r2*(Best_pos(j)-X(index2(i),j))+(1-r2)*(X_centroid(j)-X(index2(i),j)));
        end
    end
    end
    
    % Check if solutions go outside the search spaceand bring them back
    for i=1:size(X,1)
        for j=1:dim
            if X(i,j)>ub(j)
                X(i,j)=ub(j);
            end
            if X(i,j)<lb(j)
                X(i,j)=lb(j);
            end
        end
        
        % Calculate the objective values
        Objective_values(1,i)=fobj(X(i,:));
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Best_score
            Best_pos=X(i,:);
            Best_score=Objective_values(1,i);
        end
    end
   
    %% Update the elite pool
    [~,idx1]=sort(Objective_values);
    second_best=X(idx1(2),:);
    third_best=X(idx1(3),:);
    sum1=0;
    for i=1:N1
        sum1=sum1+X(idx1(i),:);
    end
    half_best_mean=sum1/N1;
    Elite_pool(1,:)=Best_pos;
    Elite_pool(2,:)=second_best;
    Elite_pool(3,:)=third_best;
    Elite_pool(4,:)=half_best_mean;

    Convergence_curve(l)=Best_score;
    l=l+1;
end
