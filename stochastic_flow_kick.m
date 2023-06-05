%%%%%%%%%%%%%%%% stochastic_flow_kick.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Summary: the following script is used to do Monte Carlo simulation of a
%   stochastic flow-kick system derived from an ODE system defined below.
%   The kicks are uniformly sampled from a given interval and the
%   inter-kick times are either uniformly sampled from a given interval or
%   are given according to a Poisson process with a specified rate. Each
%   trial of the simulation is independent and simulated using the trial()
%   function defined below. The main quantities of interest for each trial
%   are whether an excursion occured (boolean) and the time at which an
%   excursion occured (in the interval [0,infinity]).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Stochastic parameters %%%%%%%%%%%%%%
%%%% the following parameters can be modified freely
% kicks are chosen uniformly from the interval bounded by min/max_kick_size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_kick_size = 12000;
max_kick_size = 13000;
% the boolean variable exp_uniform is used to set
% whether exponential (1) or uniform (0) inter-kick times are used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_uniform = 1;
% kicks arrive according to Poisson process with rate lambda; the
% inter-kick times are exponentially distributed with parameter lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_interkick_time = 150;
lambda = 1/mean_interkick_time;
% or kicks are uniform on the interval bounded by min/max_kick_size
min_interkick_time = 7;
max_interkick_time = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Monte Carlo simulation %%%%%%%%
number_of_trials = 100;
excursion_indicator = zeros(number_of_trials,1);
excursion_time = zeros(number_of_trials,1);

tic;
parfor k = 1:number_of_trials
    [excursion_indicator(k), excursion_time(k)] = trial(exp_uniform, ...
        min_kick_size, max_kick_size, lambda, ...
        min_interkick_time, max_interkick_time);
end
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Save data %%%%%%%%%%%%%%%%%%%%%%
if exp_uniform == 1
    save(strcat("excursion_poisson", datestr(now, 'yyyy_mm_dd_HH_MM_SS'), ...
        ".mat"), "excursion_indicator", "excursion_time", "mean_interkick_time")
else
    save(strcat("excursion_uniform", datestr(now, 'yyyy_mm_dd_HH_MM_SS'), ...
        ".mat"), "excursion_indicator", "excursion_time", ...
        "min_interkick_time", "max_interkick_time")
end

%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%
histogram(excursion_time(~isnan(excursion_time)), 9, 'Normalization','PDF')

%%
%%%%%%%%%%%%%%%% trial function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Summary: the following function runs a single trial of a flow-kick
%   simulation using stochastic kicks that arrive according to a Poisson
%   process with rate lambda or are uniformly distributed on a given range;
%   the magnitude of kicks is uniformly distributed on another given range.
%   Of interest is whether a so-called excursion occurs
%   (meaning the viral load variable exceeds a predetermined threshold)
%   and, if so, the time at which such an excursion occurs.
%   Output: a boolean (0 or 1) variable indicating whether an excursion
%   occurred and a time at which the excursion occurred. When no excursion
%   occurs, the time is output as NaN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [excursion_indicator, excursion_time] = trial(exp_uniform, ...
        min_kick_size, max_kick_size, lambda, ...
        min_interkick_time, max_interkick_time)
%%%%%%%%%%%%%%%% System parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   the following constants are coefficients in the ODE model presented
%   in the paper "An immuno-epidemiological model for transient immune
%   protection: A case study for viral respi- ratory infections"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 0.35; c = 20; mu = 0.2; Beta = 5e-7; V_m = 10;
g = 0.8; C_t = 7e7; Beta_prime = 2e-5;
delta = 3; kappa = 3;
q = 1e-7; d = 2;
m_1 = 1e-4; m_2 = 0.01;
m_3 = 12000; r = 0.2; mu_prime = 0.04;
max_time = 600; % number of days to flow

%%%%%%%%%%%%%%%% ODE model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%V = Y(1), T = Y(2), I = Y(3), F = Y(4), B = Y(5), A = Y(6)
f=@(t,Y) [p*Y(3)-c*Y(1)-mu*Y(1)*Y(6)-Beta*Y(1)*Y(2)*(Y(1)/(V_m+Y(1)));...
    g*Y(2)*(1-(Y(2)+Y(3))/C_t)-Beta_prime*Y(1)*Y(2)*(Y(1)/(V_m+Y(1)));...
    Beta_prime*Y(1)*Y(2)*(Y(1)/(V_m+Y(1)))-delta*Y(3)-kappa*Y(3)*Y(4);...
    q*Y(3)-d*Y(4);...
    m_1*Y(1)*(1-Y(5))-m_2*Y(5);...
    m_3*Y(5)-r*Y(6)-mu_prime*Y(1)*Y(6)] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%
V_0 = 15000; T_0 = C_t; I_0 = 0; F_0 = 0; B_0 = 0; A_0 = 0;
IC = [V_0; T_0; I_0; F_0; B_0; A_0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Stochastic sampling %%%%%%%%%%%%%%%%%%%%%%
%%% the number of kicks is stored in the variable N
%%% inter-kick times are stored in the vector tau
if exp_uniform == 1
    N = poissrnd(max_time*lambda); % random number of kicks
    %%% kick times
    U = max_time*rand(N,1); % sample N points from interval [0,max_time]
    sorted_U = [0; sort(U); max_time]; % sort the points to get kick times
    tau = [diff(sorted_U)];
else
    N = floor(max_time/min_interkick_time);
    tau = min_interkick_time + (max_interkick_time-min_interkick_time)*rand(N,1);
    N = find(cumsum(tau) >= max_time,1) - 1; % number of kicks so that sum
    % of inter-kick times does not exceed max_time
    tau = [tau(1:N); max_time - sum(tau(1:N))];
end
%%% kick sizes
% note 0 is added as a placeholder at the end in order to match the sizes
% of the tau and kicks vectors, but 0 is never used
kicks = [min_kick_size + (max_kick_size-min_kick_size)*rand(N,1); 0];

%%%%%%%%%%%%%%% Flow kick simulation %%%%%%%%%%%%%%%%%%%%%%%%
%%% storage for output of ODE45 through all flow-kick steps
tall = [];
Yall = [];

if N>0
    for j = 1:(N+1)
        [ts, Ys] = ode45(f,[0,tau(j)], IC); % flow
        k = [kicks(j), 0, 0, 0, 0, 0]; % kicks are only applied to V
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Check trial for excursion and time of excursion %%%%%%%%
excursion_threshold = 2.5e5; % viral load considered to cause illness
minimum_time_index = find(tall >= 5,1); % first position where time is >=5

%check if viral load is ever greater than or equal to excursion_threshold
if max(Yall(minimum_time_index:end,1)) >= excursion_threshold
    excursion_time_index = find(Yall(minimum_time_index:end,1) >= ...
        excursion_threshold,1) + minimum_time_index;
    excursion_indicator = 1; %1 for excursion
    excursion_time = tall(excursion_time_index);
else
    excursion_indicator = 0; %0 for no excursion
    excursion_time = NaN;
end
end