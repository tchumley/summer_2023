load("~/Desktop/rectangle_position_vertical.mat files/rectangle_position_vertical.mat")
%plot(mean_time_to_excursion vs. lambda value) and error bar
%Then we will use fit_logistic to analyze data, find logistic regression
%expression
%presentation include: background, results(two plots, QQ-plot), logic behind code(SLLN, Pois,
%Exp, Unif)

%find probability_of_excursion and find mean_time_to_excursion
% mean(data(1).times)
% mean(data(1).results)

%write a loop to compute (32 mean_time_to_excursion) and (32
%probability_of_excursion)
mean_time_to_excursion = zeros(length(data9),1);
probability_of_excursion = zeros(length(data9), 1);
lower_bounds = zeros(length(data9),1);
std_time_to_excursion = zeros(length(data9),1);
std_probability_of_excursion = zeros(length(data9),1);
number_of_excursion = zeros(length(data9),1);

for i = 1:length(data9)
    mean_time_to_excursion(i) = mean(data9(i).times,'omitnan');
    std_time_to_excursion(i) = std(data9(i).times,'omitnan');
    probability_of_excursion(i) = mean(data9(i).indicators);
    std_probability_of_excursion(i) = std(data9(i).indicators);
    lower_bounds(i) = data9(i).lower_bounds;
    number_of_excursion(i) = mean(data9(i).indicators)*10000;
end

[lower_bounds,I] = sort(lower_bounds);
mean_time = mean_time_to_excursion(I);
probability = probability_of_excursion(I);
std_probability_of_excursion = std_probability_of_excursion(I);
std_time_to_excursion = std_time_to_excursion(I);
number_of_excursion = number_of_excursion(I);

figure(1)
plot(lower_bounds, probability, '*-')
xlabel('kick size')
ylabel('probability of excursion')


[ Qpre, p, sm, varcov] = fit_logistic(lower_bounds,probability)
p

% error bar for plot(lower bounds, probability, '*-')
upper_bound_probability = + std_probability_of_excursion./sqrt(10000)*1.96;
figure
errorbar(lower_bounds, probability, upper_bound_probability,'k')
xlabel('kick size')
ylabel('probability of excursion')
%relatviely large error bar

function [ Qpre, p, sm, varcov] = fit_logistic(t,Q)
%fit a logistic function to time series Q(t).
%   Inputs: t (time),Q (time series variable)
%   Outputs: Qpre (logistic model fit to data) and
%   p is 3 element vector containing parameters describing the logistic:
%       thalf, Qinf, and alpha
%   Q(t) = Qinf/(1 + exp(-alpha*(t-thalf)))
%       thalf is symmetric inflection point
%       Qinf is value as t --> infinity
%       alpha is time decay constant
%   sm is 3 element vector giving 1-sigma confidence limits of parameters
%       e.g., thalf = p(1) +/- sm(1)
%       simply double the values in sm to get 95% confidence limits
%   varcov is the complete 3x3 variance-covariance matrix for investigating
%       how model paramters co-vary with each other.
%       sm is sqrt(diag(varcov))
%
%   Example:
%       Qinf = 10.2; alpha = 0.33; thalf = 108.5;
%       t = 100:120;
%       Q = Qinf./(1 + exp(-alpha*(t-thalf)));
%       noise = randn(1,length(t));
%       Qd = Q+noise;
%       Qpre = fit_logistic(t,Qd);
%       figure(1)
%       clf
%       hold on
%       plot(t,Qd,'o')  % data
%       plot(t,Qpre)    % best fitting logistic
%       
% Written by James Conder, Southern Illinois University, Oct. 2010
% Cleaned up for publishing May 16, 2013
%   May 17, 2013: Allow for decreasing logistic.
%   May 23, 2013: Fix instability when using short 
%       or long absolute times (relative to alpha = 1).
%   May 28, 2013: added example in comments, fixed an introduced bug
%       from May 23 edit. 
%   Feb 12, 2014: Revisited occasional flatlining problem.
%       (Qpre goes to mean).
%       Made initial alpha more robust. Scaled to time rather than simply
%       defaulting to one (removes much of need for rescaling time).
%       Added check for flatlining. If occurs, reset seeds with larger
%       alpha and start over.
%   Jan 25, 2016: calculate confidence limits for parameters
% equations are set up to solve for an increasing logistic.
% Check if decreasing and if so, reverse time
[~,I] = sort(t);
reverse_t = false;
if sum(diff(Q(I))) < 0      % decreasing in time
    reverse_t = true;
    t = -t;
end
% stretch short or long sequences in time to stabilize alpha
tstretch = [];
if max(t)-min(t) < 1.e-4 || max(t)-min(t) > 1e5;
    tstretch = 1./(max(t) - min(t));
    t = t*tstretch;
end
% initial guesses for parameters
thalf = 0.5*(min(t) + max(t));
Qinf = max(Q);
alpha = 1./(max(t)-min(t)); alphareset = alpha;
flipQ = false;
if isrow(Q)
    flipQ = true;   % expecting a column vector. flip if row.
    t = t';
    Q = Q';
end
itermax = 1000 ;	% number of maximum iterations
epsilon = 1;
ii = 0 ;            % initialize counter
thresh = 1.e-6 ;    % threshold to stop iterating
G = zeros(length(t),3) ;    % dimensionalize partial derivative matrix
while epsilon > thresh
  ii = ii + 1 ;
  Qpre = Qinf./(1 + exp(-alpha*(t-thalf))) ;  % 'predicted' data
  if max(Qpre) - min(Qpre) == 0
      % if Qpre flatlines, "a" likely needed to be seeded higher
      % (sharper climb)
      alphareset = 2*alphareset;
      thalf = 0.5*(min(t) + max(t));
      Qinf = max(Q);
      alpha = alphareset;      
      Qpre = Qinf./(1 + exp(-alpha*(t-thalf))) ;
  end
  d = Q - Qpre ;      % data vector (predicted - observed)
  % linearized partial derivatives
  ee = min(exp(-alpha*(t-thalf)),1.e12) ;
  eee = 1./((1 + ee).^2) ;
  G(:,1) = -Qinf*alpha*(ee.*eee) ;          % dd/dthalf
  G(:,2) = 1./(1 + ee) ;                % dd/dQinf
  G(:,3) = Qinf*(t-thalf).*(ee.*eee) ;  % dd/dalpha
  
  [U,S,V] = svd(G,0);				% Singular Value Decomposition
  Sinvdiag = 1./diag(S) ;
  ising = Sinvdiag(1)./Sinvdiag < 1.e-12 ;
  Sinvdiag(ising) = 0;
  Sinv = diag(Sinvdiag);
  dm = 0.1*V*Sinv*U'*d;
							  
  % get new parameters: m = m0 + dm
  thalf = thalf + dm(1);
  Qinf = Qinf + dm(2);
  alpha = alpha + dm(3);				 
							   				  
  epsilon = norm(dm);
  if ii > itermax
	    disp('max number of iterations reached...exiting')
        disp(['normalized epsilon: ' num2str(epsilon/thresh)]) 
        epsilon = thresh ;
  end
end
Qpre = Qinf./(1 + exp(-alpha*(t-thalf))) ;  % 'predicted' data
if ~isempty(tstretch)
    thalf = thalf/tstretch;
    alpha = alpha*tstretch;
end  
if reverse_t        % decreasing logistic
    thalf = -thalf;
    alpha = -alpha;
end
if flipQ
    Qpre = Qpre';
    Q = Q';     % necessary for a posteriori variance 
end
p = [ thalf Qinf alpha ];
%%% find confidence bounds for parameters (1-sigma)
if nargout > 2
    sd = sum((Q-Qpre).^2)/(length(Q)-3);       % a posteriori variance
    varcov = sd*sd*inv(G'*G);   % model variance-covariance matrix
    sm = sqrt(diag(varcov));    % model variance
end
end


