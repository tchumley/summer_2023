%cd('Desktop')
cd("~/Desktop/rectangle_position_horizontal_8000/")
files = dir('*.mat');
times = zeros(length(files), 10000);
indicators = zeros(length(files), 10000);
mean_times = zeros(length(files), 1);
probabilities = zeros(length(files), 1);
lower_bounds = zeros(length(files), 1);
taus = zeros(length(files), 1);


for i = 1:length(files)
    load(files(i).name);
    times(i,:) = excursion_time;
    indicators(i,:) = excursion_indicator;
    mean_times(i) = mean_time_to_excursion;
    lower_bounds(i) = lower_bound;
    probabilities(i) = probability_of_excursion;
    taus(i) = which_tau;
end

% our_lambdas = 7:3:100;
% 
% for i = 1:length(our_lambdas)
%     indices = find(lambdas == our_lambdas(i))/9;
%     times(i,1:10000) = times(indices,:);
%     indicators(i,1:10000) = indicators(indices,:);
% end
% 
% for k = 1:length(our_lambdas)
%     data(k).lambda = our_lambdas(k);
%     data(k).times = times(k,:);
%     data(k).indicator = indicators(k,:);
% end

for i = 1:length(files)
    data15(i).times = times(i,:);
    data15(i).indicators = indicators(i,:);
    data15(i).mean_times = mean_times(i);
    data15(i).lower_bounds = lower_bounds(i);
    data15(i).probabilities = probabilities(i);
    data15(i).taus = taus(i);
end

save("./rectangle_position_horizontal_8000.mat", 'data15')