load("first_excursion_time_1_10000.mat")
histogram(inter_excursion_time)
histfit(inter_excursion_time,20,'exponential')
xlabel('time')
ylabel('number of trials')
pd = fitdist(inter_excursion_time,'exponential')

qqplot(inter_excursion_time,pd)
grid on
hold on
box on

% histogram(T_combine)
% histfit(T_combine,20,'exponential')
% title('histogram of the excursion time(T), mean of tau=4, max time=10000, threshold=50')
% xlabel('time')
% ylabel('number of trials')
% pd = fitdist(T_combine,'exponential')

% histogram(sigma_combine)
% histfit(sigma_combine,20,'exponential')
% title('histogram of the inter limit cycle time(sigma), mean of tau=4, max time=10000, threshold=50')
% xlabel('time')
% ylabel('number of trials')
% pd = fitdist(sigma_combine,'exponential')
%

% histogram(T);
% % p.FaceColor = "#0000FF";
% % p.EdgeColor = 'none';
% % histfit(T,20,'exponential')
% h = histfit(T,20,'exponential');
% set(h(1),'facecolor',"#0072BD"); set(h(2),'color','none')
% % title('histogram of the excursion time(T), mean of tau=4, max time=10000, threshold=50')
% xlabel('time')
% ylabel('number of trials')
% grid on
% % % pd = fitdist(T,'Normal')
% % % qqplot(T,pd)

% histogram(sigma)
% histfit(sigma,20,'exponential')
% % title('histogram of the inter limit cycle time(sigma), mean of tau=4, max time=10000, threshold=50')
% xlabel('time')
% ylabel('number of trials')
% grid on
% % pd = fitdist(sigma,'exponential')
% % qqplot(sigma,pd)