number_of_trials = 1000;
trail_number = 1:number_of_trials;
results = zeros(number_of_trials,1);
for k = 1:number_of_trials
    results(k) = trial();
end
mean(results)


function trial_result = trial()
f=@(t,Y) [Y(1)*(1-Y(1))*(Y(1)-2)-Y(2);...
    0.1*(Y(1)-0.2)];
IC = [0.2;-0.288];
max_time = 10000;
mean_interkick_time = 10;
lambda = 1/mean_interkick_time;
N = poissrnd(max_time*lambda); % random number of kicks
%%% kick times
U = max_time*rand(N,1); % sample N points from interval [0,max_time]
sorted_U = [0; sort(U); max_time]; % sort the points to get kick times
tau = [diff(sorted_U)];

min_kick_size = -.1;
max_kick_size = 0;
kicks = [min_kick_size + (max_kick_size-min_kick_size)*rand(N,1); 0];

[t,y] = ode45(f,[0,max_time],IC);

tall = [];
Yall = [];

if N>0
    for j = 1:(N+1)
        [ts, Ys] = ode45(f,[0,tau(j)], IC); % flow
        k = [0, kicks(j)]; % kicks are only applied to predators
        IC = Ys(end,:) + k; % add kick for next flow
        % store times output by ODE45
        % (shifted by previous total kick times)
        % and store ODE variables output by ODE45
        tall = [tall; ts+sum(tau(1:(j-1)))]; 
        Yall = [Yall; Ys];
    end
else % special case when no kicks occur
    [tall, Yall] = ode45(f,[0,max_time], IC);
end
% figure
% plot(Yall(:,1), Yall(:,2))
% hold on
% fplot(@(x) x.*(1-x).*(x-2),[-1,3])
% 
% figure
% plot(tall,Yall)
% 
trial_result = excursioncheck(tall, Yall);
end


function excursion_indicator = excursioncheck(tall,Yall)
ymax = max(Yall);
% ymax(1);
excursion_threshold = 1+1/sqrt(3);
if ymax(1) >= excursion_threshold
    excursion_indicator = 1; %1 for excursion
else
    excursion_indicator = 0; %0 for no excursion
end
end