% peak detector adapted from http://billauer.co.il/blog/2009/01/peakdet-matlab-octave/
f=@(t,Y) [Y(1)*(1-Y(1))*(Y(1)-2)-Y(2);...
    0.1*(Y(1)-0.2)];
%%%%%%%%% Stochastic parameters %%%%%%%%%%%%%%
IC = [0.2;-0.288];
max_time = 1000;
mean_interkick_time = 10;
lambda = 1/mean_interkick_time;
min_kick_size = -0.1;
max_kick_size = 0;
number_of_excursion = 0;
filled = 0;
threshold = 50;
number_of_trials = 1;
excursion_index = 2;
list_index = 1;
inter_peak_time = zeros(2000,1);
inter_excursion_time = zeros(2000,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%% Generate Inter-peak-time %%%%%%
for k = 1:number_of_trials   
    [tall,Yall] = trial(max_time, ... 
            min_kick_size, max_kick_size, lambda, f, IC);
    t=tall;
    x=Yall(:,1);
    [maxtab, mintab] = peakdet(x, 1, t);

    if list_index == 1
        inter_peak_time(list_index) = maxtab(1,1);
        list_index = list_index + 1;
        for j = 1:length(mintab)
                if (mod(list_index,2) == 0) && (list_index+1 <= length(inter_peak_time))
                inter_peak_time(list_index) = mintab(j,1)-maxtab(j,1);
                list_index = list_index + 1;
                    if (j+1 <= length(maxtab)) && (list_index+1 <= length(inter_peak_time))
                        inter_peak_time(list_index) = maxtab(j+1,1)-mintab(j,1);
                        list_index = list_index + 1;
                    end
                else
                    break
                end
        end
    else
        for j = 1:length(mintab)
                if (mod(list_index,2) == 0) && (list_index+1 <= length(inter_peak_time))
                inter_peak_time(list_index) = mintab(j,1)-maxtab(j,1);
                list_index = list_index + 1;
                    if (j+1 <= length(maxtab)) && (list_index+1 <= length(inter_peak_time))
                        inter_peak_time(list_index) = maxtab(j+1,1)-mintab(j,1);
                        list_index = list_index + 1;
                    end
                else
                    break
                end
        end
    end
end
inter_peak_time = nonzeros(inter_peak_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; plot(t,x);
hold on; plot(mintab(:,1), mintab(:,2), 'g*');
plot(maxtab(:,1), maxtab(:,2), 'r*');


%%% Generate Inter-excursion-time %%%
inter_excursion_time(1) = inter_peak_time(1);
for j = 2:length(inter_peak_time)    
    if (inter_peak_time(j) <= threshold)
        while (inter_peak_time(j) <= threshold) && (excursion_index <= length(inter_excursion_time)) 
            inter_excursion_time(excursion_index) = inter_excursion_time(excursion_index) + inter_peak_time(j);
            if j+1 <= length(inter_peak_time)
                j=j+1;
            else
                break
            end
        end
        excursion_index = j;
    else
        inter_excursion_time(excursion_index) = inter_peak_time(excursion_index);
        excursion_index = excursion_index + 1;
    end
end
        
inter_excursion_time = nonzeros(inter_excursion_time);
number_of_excursion = length(inter_excursion_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%% Generate T and sigma %%%%%%%%
sigma = zeros(fix(number_of_excursion/2),1);
T = zeros(fix(number_of_excursion/2),1);
s_index = 1;
T_index = 1;
        
for l = 1:length(inter_excursion_time)
    if mod(l,2) ~= 0       
        sigma(s_index,1) = inter_excursion_time(l);
        s_index = s_index + 1;
    else
        T(T_index,1) = inter_excursion_time(l);
        T_index = T_index + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save(strcat("inter_excursion_time", datestr(datetime), ".mat"),"number_of_excursion","inter_peak_time", "inter_excursion_time", "sigma", "T")




%%%%%%%%%%%%%%%% peakdet function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05
% This function is released to the public domain; Any use is allowed.

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
  x = (1:length(v))';
else
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end

if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(v)
  this = v(i);
  if this > mx, mx = this; mxpos = x(i); end
  if this < mn, mn = this; mnpos = x(i); end

  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end
  else
    if this > mn+delta
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end
end
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
function [tall,Yall] = trial(max_time, ... 
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
