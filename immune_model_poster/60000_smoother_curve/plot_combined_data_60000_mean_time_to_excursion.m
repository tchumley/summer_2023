load("~/Desktop/60000/data_combined60000.mat")
%plot(mean_time_to_excursion vs. lambda value) and error bar
%plot(mean_excursion_results vs. lambda value) and error bar
%Then we will use QQ-plot to analyze data
%presentation include: background, results(two plots, QQ-plot), logic behind code(SLLN, Pois,
%Exp, Unif)

%find probability_of_excursion and find mean_time_to_excursion
% mean(data(1).times)
% mean(data(1).results)

%write a loop to compute (32 mean_time_to_excursion) and (32
%probability_of_excursion)
mean_time_to_excursion = zeros(length(data13),1);
probability_of_excursion = zeros(length(data13), 1);
lambdas = zeros(length(data13),1);
std_time_to_excursion = zeros(length(data13),1);
std_probability_of_excursion = zeros(length(data13),1);

for i = 1:length(data13)
    mean_time_to_excursion(i) = mean(data13(i).times,'omitnan');
    std_time_to_excursion(i) = std(data13(i).times,'omitnan');
    probability_of_excursion(i) = mean(data13(i).indicators);
    std_probability_of_excursion(i) = std(data13(i).indicators);

    lambdas(i) = data13(i).lambdas;
end

[lambdas,I] = sort(lambdas);
mean_time = mean_time_to_excursion(I);
probability = probability_of_excursion(I);
std_probability_of_excursion = std_probability_of_excursion(I);
std_time_to_excursion = std_time_to_excursion(I);

% figure(1)
% plot(1./lambdas, probability, '*-')
% xlabel('tau')
% ylabel('probability of excursion')
% 
% figure(2)
% plot(lambdas, probability, '*-')
% xlabel('lambda (average number of kicks in 1 days)')
% ylabel('probability of excursion')
% 
% figure(3)
% plot(1./lambdas, mean_time, '*-')
% xlabel('mean of tau')
% ylabel('mean time to excursion')
% 
% 
% figure(4)
% plot(lambdas, mean_time, '*-')
% xlabel('lambda (average number of kicks in 1 days)')
% ylabel('mean time to excursion')



% error bar for plot(1./lambdas, probability, '*-')
upper_bound_probability_lambdas = + std_probability_of_excursion./sqrt(60000)*1.96;
figure(1)
errorbar(1./lambdas, probability, upper_bound_probability_lambdas,'k')
xlabel('mean of \tau')
ylabel('probability of excursion')
yticks([0.95 0.96 0.97 0.98 0.99 1])
ax = gca; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on
set(gca,'LooseInset',get(gca,'TightInset'));
%title('probability of excursion vs. mean of tau, n = 60000','FontSize',16)

% %relatviely large error bar
% 
% %error bar for plot(lambdas, probability, '*-')
% upper_bound_probability = + std_probability_of_excursion./sqrt(60000)*1.96;
% figure
% errorbar(lambdas, probability, upper_bound_probability)
% xlabel('lambda (average number of kicks in 1 days)')
% ylabel('probability of excursion')

%the errorbar is very tiny

%error bar for plot(1./lambdas, mean_time, '*-')
upper_bound_mean_time_lambdas = + std_time_to_excursion./sqrt(probability_of_excursion*60000)*1.96;
figure(2)
errorbar(1./lambdas, mean_time, upper_bound_mean_time_lambdas,'k')
xlabel('mean of tau')
ylabel('mean time to excursion')
xlim([0 100])
ax = gca; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on
set(gca,'LooseInset',get(gca,'TightInset'));
title('mean time to excursion vs. mean of tau, n = 60000','FontSize',16)
%the errorbar is very tiny, so we will not analyze it

% %error bar for plot(lambdas, mean_time, '*-')
% upper_bound_mean_time = + std_time_to_excursion./sqrt(60000)*1.96;
% figure
% errorbar(lambdas, mean_time, upper_bound_mean_time)
% xlabel('lambda (average number of kicks in 1 days)')
% ylabel('mean time to excursion')
% %the errorbar is very tiny, so we will not analyze it

eps_save(1,'60000probability of excursion')

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
