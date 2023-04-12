clear all
clc
close all
SearchAgents_no=30; % Number of search agents
% Function_name='F12'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration=1000; % Maximum number of iterations  
Function_name=1; 
dim=10; %Dimension,[2,10,30,50,100]
lb=-100;%lower boound
ub=100;%upper bound
fobj = @(x) cec17_func(x',Function_name);

Max_test=20;
for i=1:Max_test
    disp(['The',num2str(i),' -th experiment']);
    [Best_pos(i,:),Best_score(i),SAO_curve(i,:)]=SAO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化
end

figure
semilogy(mean(SAO_curve),'color','k','linewidth',2.0,'Marker','o','MarkerIndices',1:100:length(mean(SAO_curve)))
title('Convergence curve of C2017_{1}')
xlabel('Iteration');
ylabel('Fitness');
axis tight
grid off
box on 
legend('SAO')

disp('-------------------------------------------------')
display(['The best fitness value obtained by SAO over 20 independent runs : ', num2str(min(Best_score))]);
display(['The mean fitness value obtained by SAO over 20 independent runs : ', num2str(mean(Best_score))]);
display(['The worst fitness value obtained by SAO over 20 independent runs : ', num2str(max(Best_score))]);
display(['The standard deviation obtained by SAO over 20 independent runs : ', num2str(std(Best_score))]);

