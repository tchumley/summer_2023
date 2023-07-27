load("~/Desktop/mean_tau_vs_probability.mat")
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
probability_of_excursion = zeros(length(data10), 1);
mean_of_tau = zeros(length(data10),1);
std_probability_of_excursion = zeros(length(data10),1);

for i = 1:length(data10)
    probability_of_excursion(i) = mean(data10(i).indicators);
    std_probability_of_excursion(i) = std(data10(i).indicators);
    mean_of_tau(i) = data10(i).mean_tau;
end

[mean_of_tau,I] = sort(mean_of_tau);
probability = probability_of_excursion(I);
std_probability_of_excursion = std_probability_of_excursion(I);



% error bar for plot(1./lambdas, probability, '*-')
upper_bound_probability_tau = + std_probability_of_excursion./sqrt(20000)*1.96;
figure
errorbar(mean_of_tau, probability, upper_bound_probability_tau,'k')
xlabel('mean of tau')
ylabel('probability of excursion')
grid on
% title('max time = 1000')
%relatviely large error bar

% %error bar for plot(lambdas, probability, '*-')
% upper_bound_probability = + std_probability_of_excursion./sqrt(20000)*1.96;
% figure
% errorbar(lambdas, probability, upper_bound_probability)
% xlabel('lambda (average number of kicks in 1 days)')
% ylabel('probability of excursion')
% 
% %the errorbar is very tiny
% 
% %error bar for plot(1./lambdas, mean_time, '*-')
% upper_bound_mean_time_lambdas = + std_time_to_excursion./sqrt(100000)*1.96;
% figure
% errorbar(1./lambdas, mean_time, upper_bound_mean_time_lambdas,'k')
% xlabel('mean of tau')
% ylabel('mean time to excursion')
% %the errorbar is very tiny, so we will not analyze it
% 
% %error bar for plot(lambdas, mean_time, '*-')
% upper_bound_mean_time = + std_time_to_excursion./sqrt(100000)*1.96;
% figure
% errorbar(lambdas, mean_time, upper_bound_mean_time)
% xlabel('lambda (average number of kicks in 1 days)')
% ylabel('mean time to excursion')
% %the errorbar is very tiny, so we will not analyze it
