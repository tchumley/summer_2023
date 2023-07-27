combined = importdata("rectangle_data_combine.mat");

% mean_time_to_excursion = zeros(length(combined),1);
% probability_of_excursion = zeros(length(combined), 1);
% lower_bounds = zeros(length(combined),1);
% taus = zeros(length(combined),1);
% std_time_to_excursion = zeros(length(combined),1);
% std_probability_of_excursion = zeros(length(combined),1);
% number_of_excursion = zeros(length(combined),1);

%%
% for i = 1:length(combined)
%     mean_time_to_excursion(i) = mean(combined(i).times,'omitnan');
%     std_time_to_excursion(i) = std(combined(i).times,'omitnan');
%     probability_of_excursion(i) = mean(combined(i).indicators);
%     std_probability_of_excursion(i) = std(combined(i).indicators);
%     taus(i) = combined(i).taus;
%     lower_bounds(i) = combined(i).lower_bounds;
%     number_of_excursion(i) = mean(combined(i).indicators)*10000;
% end

% [taus,I] = sort(taus);
% mean_time = mean_time_to_excursion(I);
% probability = probability_of_excursion(I);
% std_probability_of_excursion = std_probability_of_excursion(I);
% std_time_to_excursion = std_time_to_excursion(I);
% number_of_excursion = number_of_excursion(I);


%plot3(lower_bounds,taus,probability_of_excursion )

probabilities = zeros(11, 9);
kvalue = [1000:1000:11000];
tvalue = [4:1:12];

for j = 1:length(combined)
    tindex = find(combined(j).taus == tvalue);
    bound_index = find(combined(j).lower_bounds == kvalue);
    probabilities(bound_index,tindex) = combined(j).probabilities;

end

%%
figure(1)
[t,k] = meshgrid(tvalue,kvalue);

surf(t,k,probabilities)
xlabel('mean flow time')
ylabel('kick size')
zlabel('probability of excursion')

%%
figure(2)
contourf(tvalue,kvalue,probabilities)
xlabel('mean flow time')
ylabel('kick size')
colorbar
clim([0 1])
