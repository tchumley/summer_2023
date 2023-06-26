alpha = 0.01; beta = 0.02;
f=@(t,Y) [(1-alpha*Y(2))*Y(1);...
    (-1+beta*Y(1))*Y(2)];

IC = [10;10];

max_time = 1000;
mean_interkick_time = 10;
lambda = 1/mean_interkick_time;
N = poissrnd(max_time*lambda); % random number of kicks
%%% kick times
U = max_time*rand(N,1); % sample N points from interval [0,max_time]
sorted_U = [0; sort(U); max_time]; % sort the points to get kick times
tau = [diff(sorted_U)];

min_kick_size = 50;
max_kick_size = 100;
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

figure
%plot(tall,Yall)
%hold on
%plot(t,y)
plot(y(:,1),y(:,2), '*')
hold on
%figure
plot(Yall(:,1), Yall(:,2))