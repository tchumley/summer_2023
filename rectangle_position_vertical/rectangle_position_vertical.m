% %load('final_results_uniform.mat')
% %this gives a 1x1000 vector of 0s and 1s
% A = cumsum(results);
% probability_of_trails = cumsum(results)./(1:1000);
% plot(1:1000, probability_of_trails);
% %first define mu, sigma, n, and the range for the 1000th trail.
% mu = mean(results);
% sigma = std(results);
% n = 1000;
% a = mu + sigma/sqrt(n)*1.96;
% b = mu - sigma/sqrt(n)*1.96;
% data = [mu, sigma];
% save('mu_sigma.mat', 'data')
% %second construct errorbar.
% x = floor(linspace(1,1000,21));
% y = probability_of_trails(floor(linspace(1,1000,21)));
% upper_bound = + sigma./sqrt(x)*1.96;
% lower_bound = - sigma./sqrt(x)*1.96;
% errorbar(x(5:21),y(5:21),upper_bound(5:21),lower_bound(5:21))
%mean time to excursion
%sum = zeros(20,1);
%% you don't need this part anymore - you are already did it in trail_2_20_trial.m file
%V = Y(1), T = Y(2), I = Y(3), F = Y(4), B = Y(5), A = Y(6)
rng('shuffle');
number_of_trials = 10000;
trail_number = 1:number_of_trials;
excursion_indicator = zeros(number_of_trials, 1);
excursion_time = zeros(number_of_trials, 1);

tic; 
parfor k = 1:number_of_trials %number_of_trials
    %[results(k,1),results(k,2)] = trial();
    [excursion_indicator(k),excursion_time(k)] = trial();
end
toc;
list_of_ones = find(excursion_indicator(:) == 1);
mean_time_to_excursion = sum(excursion_time(list_of_ones)./length(list_of_ones));
probability_of_excursion = mean(excursion_indicator);
lower_bound = 13000;
 
save(strcat("different_rectangle_position_uniform", datestr(datetime), ".mat"), "lower_bound", "excursion_indicator", "excursion_time","mean_time_to_excursion", "probability_of_excursion")

%% 
%get rid of NaN
%results(1:20,1) = randi([0,1],[1,20])
%results(1:20,2) = 18*rand(20,1)
%find(results(:,1) == 1)
%list_of_ones = find(excursion_result(:) == 1);
%results(list_of_ones,:)
%length(list_of_ones)
%mean_time_to_excursion = sum(excursion_time(list_of_ones)./length(list_of_ones));
%mean(results(:,2))
%histogram
%histogram(excursion_time(list_of_ones), 9)

%A = cumsum(results)
%probability_of_trails = cumsum(results)./trail_number;
%plot(trail_number, probability_of_trails);

%mean(results)
%save('results5.mat','results')


function [excursion_indicator,excursion_time] = trial()
f=@(t,Y) [0.35*Y(3)-20*Y(1)-0.2*Y(1)*Y(6)-5*10^(-7)*Y(1)*Y(2)*(Y(1)/(10+Y(1)));...
    0.8*Y(2)*(1-(Y(2)+Y(3))/(7*10^7))-2*10^(-5)*Y(1)*Y(2)*(Y(1)/(10+Y(1)));...
    2*10^(-5)*Y(1)*Y(2)*(Y(1)/(10+Y(1)))-3*Y(3)-3*Y(3)*Y(4);...
    1*10^(-7)*Y(3)-2*Y(4);...
    1*10^(-4)*Y(1)*(1-Y(5))-0.01*Y(5);...
    12000*Y(5)-0.2*Y(6)-0.04*Y(1)*Y(6)] ;

IC = [15000; 7*10^7; 0; 0; 0; 0];
k_V = 0;
k = [k_V 0 0 0 0 0];

tau = 7+1*rand(85,1);
steps = find(cumsum(tau) >= 600,1)-1;
tau = tau(1:steps);
a = 0;
b = 1000;
tall = [];
Yall = [];
kall = [];

for c = 1:steps
    [ts, Ys] = ode45(f,[0,tau(c)], IC);
    k = [13000 + (a+(b-a)*rand), 0, 0, 0, 0, 0];
    IC = Ys(end,:) + k;
    tall = [tall; ts+sum(tau(1:(c-1)))];
    Yall = [Yall; Ys];
    kall = [kall; k(1)];
    %excursioncheck(tall, Yall);
end
%plot(tall,Yall(:,1))
[excursion_indicator,excursion_time] = excursioncheck(tall, Yall);
end



function [excursion_indicator,excursion_time] = excursioncheck(ts,Ys)
minimum_time = find(ts >= 5,1); %first position where time is >=5
if max(Ys(minimum_time:end,1)) >= 2.5*10^5 %max of virsus in positions minimum_time to end, checking if that is >=2*10^5
    %display('excursion')
    index_of_excursion = find(Ys(minimum_time:end,1) >= 2.5*10^5,1)+minimum_time;
    Ys(index_of_excursion,1);
    ts(index_of_excursion);
    excursion_indicator = 1; %1 for excursion
    excursion_time = ts(index_of_excursion);
else
    excursion_indicator = 0;
    excursion_time = NaN;
    %display ('no excursion')
end
end
