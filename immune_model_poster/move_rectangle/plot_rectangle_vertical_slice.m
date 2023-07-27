probabilities = importdata("probability_of_excursion_data_matrix.mat");
kvalue = [1000:1000:11000];
tvalue = [4:1:12];
[t,k] = meshgrid(tvalue,kvalue);

% create nine probability_matrix that record probability of excursion based on kick sizes for each tau value
% plot each probability_matrix to get a slice
% hold on to get the whole graph
figure(1)
for i = 1:9
    A=zeros(2,1)+tvalue(i);
    probability_matrix = zeros(11, 2);
    probability_matrix(:,2) = probabilities(:,i);
    surf(A,kvalue,probability_matrix,'EdgeColor','white')
    hold on
end

xlabel('mean flow time')
ylabel('kick size')
zlabel('probability of excursion')

grid on
ylim([1000 11000])
zticks([0,0.2,0.4,0.6,0.8,1.0])
ax = gca; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.ZAxis.FontSize = 16;
set(gca,'LooseInset',get(gca,'TightInset'));
%title('kick size vs. probability of excursion, different \tau slices','FontSize',16)
hold off

eps_save(1,'rectangle_3D')

function y = eps_save(fig_number,filename)
figure(fig_number)
set(gcf,'PaperUnits','inches');
oldsizes = get(gcf,'PaperPosition');
% This returns [x y width height]
newwidth = 3.2;
newheight = oldsizes(4)/oldsizes(3)*newwidth;
set(gcf,'PaperPosition',[0 0 newwidth newheight]);
print('-opengl',filename,'-depsc','-r300')
end
