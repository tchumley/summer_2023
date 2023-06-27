load("our_data.mat")

mean_time_to_excursion = zeros(length(data),1);
median_time_to_excursion = zeros(length(data), 1);
lambdas = zeros(length(data),1);
mode_time_to_excursion = zeros(length(data),1);

for i = 1:length(data)
    mean_time_to_excursion(i) = mean(data(i).times,'omitnan');
    median_time_to_excursion(i) = median(data(i).times,'omitnan');
    mode_time_to_excursion(i) = mode(data(i).times);
    lambdas(i) = data(i).lambda;
end
[lambdas,I] = sort(lambdas);
mean_time = mean_time_to_excursion(I);
median = median_time_to_excursion(I);
mode = mode_time_to_excursion(I);

%when try to plot three curves onto one figure, it does not work
figure(1)
hold on
a1 = plot(600./lambdas, median, '*-')
% ylabel('median time to excursion')
xlabel('mean of tau')
ylabel('time to excursion')
M1 = "median"

a2 = plot(600./lambdas, mode, '*-')
% xlabel('mean of tau')
% ylabel('mode time to excursion')
M2 = "mode"

a3 = plot(600./lambdas, mean_time, '*-')
% xlabel('mean of tau')
% ylabel('mean time to excursion')
M3 = "mean"
legend([a1,a2,a3], [M1, M2, M3])