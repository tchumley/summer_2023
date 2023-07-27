load("data_combined60000.mat")

mean_time_to_excursion = zeros(length(data13),1);
median_time_to_excursion = zeros(length(data13), 1);
lambdas = zeros(length(data13),1);
mode_time_to_excursion = zeros(length(data13),1);

for i = 1:length(data13)
    mean_time_to_excursion(i) = mean(data13(i).times,'omitnan');
    median_time_to_excursion(i) = median(data13(i).times,'omitnan');
    mode_time_to_excursion(i) = mode(data13(i).times);
    lambdas(i) = data13(i).lambdas;
end
[lambdas,I] = sort(lambdas);
mean_time = mean_time_to_excursion(I);
median = median_time_to_excursion(I);
mode = mode_time_to_excursion(I);

%when try to plot three curves onto one figure, it does not work
figure(1)
hold on
a1 = plot(1./lambdas, median, '*-')
% ylabel('median time to excursion')
xlabel('mean of \tau')
ylabel('time to excursion')
M1 = "median"

a2 = plot(1./lambdas, mode, '*-')
% xlabel('mean of tau')
% ylabel('mode time to excursion')
M2 = "mode"

a3 = plot(1./lambdas, mean_time, '*-')
% xlabel('mean of tau')
% ylabel('mean time to excursion')
M3 = "mean"
legend([a1,a2,a3], [M1, M2, M3])
xlim([10 100])
ax = gca; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
grid on
set(gca,'LooseInset',get(gca,'TightInset'));
%title('mean \tau vs. mean, median, mode of excursion','FontSize',16)
fontsize(legend,16,"points");
hold off

eps_save(1,'mean_median_mode')

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
