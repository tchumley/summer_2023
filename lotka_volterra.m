alpha = 0.01; beta = 0.02;
f=@(t,Y) [(1-alpha*Y(2))*Y(1);...
    (-1+beta*Y(1))*Y(2)];

IC = [10;10];
max_time = 50;
[t,y] = ode45(f,[0,max_time],IC);
figure
plot(y(:,1),y(:,2))
figure
plot(t,y)