% peak detector adapted from http://billauer.co.il/blog/2009/01/peakdet-matlab-octave/
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

% t=tall;
% x=Yall(:,1);
number_of_trials = 10000;
first_excursion_time = zeros(number_of_trials,1);
excursion_indicator = zeros(number_of_trials,1);

for j = 1:length(mean_interkick_time)
    for k = 1:number_of_trials   
        [tall,Yall,excursion_indicator(k)] = trial(max_time, ... 
            min_kick_size, max_kick_size, lambda, f, IC);
        if excursion_indicator(k) == 1
            t=tall;
            x=Yall(:,1);
            [maxtab, mintab] = peakdet(x, 1, t);
%         try
%         if
            first_excursion_time(k,1) = maxtab(1,1);
%         catch
%             figure; plot(t,x);
%             ME = MException('out of bound');
%             throw(ME)
%             
%         end
        else
            first_excursion_time(k,1) = NaN;
        end
%         figure; plot(t,x);
%         hold on; 
%         plot(maxtab(:,1), maxtab(:,2), 'r*');
    end
save(strcat("inter_excursion_time", datestr(datetime), ".mat"),"first_excursion_time","max_time","mean_interkick_time","excursion_indicator")
end
inter_excursion_time(1) = maxtab(1,1);


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
      if i>1
          break
      end
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

function [tall,Yall,excursion_indicator] = trial(max_time, ... 
            min_kick_size, max_kick_size, lambda, f, IC);
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