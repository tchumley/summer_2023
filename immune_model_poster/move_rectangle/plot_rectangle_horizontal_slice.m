probabilities = importdata("probability_of_excursion_data_matrix.mat");
kvalue = [1000:1000:11000];
tvalue = [4:1:12];
[t,k] = meshgrid(tvalue,kvalue);

% create eleven probability_matrix that record probability of excursion based on tau values for each kick size
% plot each probability_matrix to get a slice
% hold on to get the whole graph
for i = 1:11
    A=zeros(1,2)+kvalue(i);
    probability_matrix = zeros(9, 2);
    probability_matrix(:,2) = probabilities(i,:);
    surf(A,tvalue,probability_matrix,'EdgeColor','white')
    hold on
end

xlabel('kick size')
ylabel('mean flow time')
zlabel('probability of excursion')

grid on
xlim([1000 11000])
ax = gca; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.ZAxis.FontSize = 16;
set(gca,'LooseInset',get(gca,'TightInset'));
title('mean flow time vs. probability of excursion, different k slices','FontSize',16)

hold off
