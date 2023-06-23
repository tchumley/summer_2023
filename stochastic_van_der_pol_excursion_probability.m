% alpha = 0.01; beta = 0.02;
f=@(t,Y) [Y(1)*(1-Y(1))*(Y(1)-2)-Y(2);...
    0.1*(Y(1)-0.2)];

%%%%%%%%% Stochastic parameters %%%%%%%%%%%%%%
IC = [0.2;-0.288];
max_time = 10000;
mean_interkick_time = 10;
lambda = 1/mean_interkick_time;
min_kick_size = -0.1;
max_kick_size = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Monte Carlo simulation %%%%%%%%
number_of_trials = 10000;
excursion_indicator = zeros(number_of_trials,1);
% excursion_time = zeros(number_of_trials,1);

for j = 1:length(mean_interkick_time)
    for k = 1:number_of_trials    
    [excursion_indicator(k)] = trial(max_time, ... 
            min_kick_size, max_kick_size, lambda, f, IC);
           
    end
%%%calculate the probability of excursion
probability_excursion = sum(excursion_indicator)/number_of_trials;
save(strcat("excursion_result_stochastic", datestr(datetime), ".mat"),"excursion_indicator", "probability_excursion")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% trial function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Summary: the following function runs a single trial of a stochastic 
%   flow-kick simulation using stochastic kicks that arrive according to
%   a Poisson process with rate lambda;
%   Of interest is whether a so-called excursion occurs
%   (meaning the viral load variable exceeds a predetermined threshold)
%   and, if so, the time at which such an excursion occurs.
%   Output: a boolean (0 or 1) variable indicating whether an excursion
%   occurred and a time at which the excursion occurred. When no excursion
%   occurs, the time is output as NaN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [excursion_indicator] = trial(max_time, ... 
            min_kick_size, max_kick_size, lambda, f, IC)
    N = poissrnd(max_time*lambda); % random number of kicks
    %%% kick times
    U = max_time*rand(N,1); % sample N points from interval [0,max_time]
    sorted_U = [0; sort(U); max_time]; % sort the points to get kick times
    tau = [diff(sorted_U)];
    kicks = [min_kick_size + (max_kick_size-min_kick_size)*rand(N,1); 0];
            
    [t,y] = ode45(f,[0,max_time],IC);
            
    tall = [];
    Yall = [];
        
    if N>0
        for j = 1:(N+1)
            [ts, Ys] = ode45(f,[0,tau(j)], IC); % flow
            kick = [0, kicks(j)]; % kicks are only applied to predators
            IC = Ys(end,:) + kick; % add kick for next flow
            % store times output by ODE45
            % (shifted by previous total kick times)
            % and store ODE variables output by ODE45
            tall = [tall; ts+sum(tau(1:(j-1)))]; 
            Yall = [Yall; Ys];
        end
        else % special case when no kicks occur
            [tall, Yall] = ode45(f,[0,max_time], IC);
    end
        %%%check if viral load is ever greater than or equal to excursion_threshold
        excursion_threshold = 1+1/sqrt(3); % viral load considered to cause illness
        if max(Yall(1:end,1)) >= excursion_threshold;
            excursion_indicator = 1; %1 for excursion
            %excursion_time(k) = tall(excursion_time_index);
        else
            excursion_indicator = 0; %0 for no excursion
            %excursion_time(k) = NaN;
        end
end


% figure
% plot(tall,Yall)
%hold on
%plot(t,y)
%plot(y(:,1),y(:,2), '*')
% figure
% fplot(@(x) x.*(1-x).*(x-2),[-1,3])
% hold on
% plot(Yall(:,1), Yall(:,2))
