load('uniform10000.mat')
A = transpose(excursion_result)
probability_of_excursion = cumsum(A)./(1:10000);
figure(1)
plot(1:10000, probability_of_excursion);
xlabel('number of trial')
ylabel('probability of excursion')
ax = gca; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on
set(gca,'LooseInset',get(gca,'TightInset'));
%title('10000 trials','FontSize',16)

% ylim([0.55,0.8])
% mu = mean(excursion_result)
% sigma = std(excursion_result)
% n = 10000
% a = mu + sigma/sqrt(n)*1.96
% b = mu - sigma/sqrt(n)*1.96
% data = [mu, sigma];
% % save('mu_sigma.mat', 'data')
% %second construct errorbar.
% x = floor(linspace(1,10000,21));
% y = probability_of_excursion(floor(linspace(1,10000,21)))
% upper_bound = + sigma./sqrt(x)*1.96
% lower_bound = - sigma./sqrt(x)*1.96
% %% 
% errorbar(x(5:21),y(5:21),upper_bound(5:21),lower_bound(5:21))

% xlabel('number of trial')
% ylabel('probability of excursion')
% ax = gca; 
% ax.XAxis.FontSize = 16;
% ax.YAxis.FontSize = 16;
% grid on
% set(gca,'LooseInset',get(gca,'TightInset'));
% title('10000 trials','FontSize',16)

eps_save(1,'converge')

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
