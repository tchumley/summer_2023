f=@(t,Y) [Y(1)*(1-Y(1))*(Y(1)-2)-Y(2);...
    0.1*(Y(1)-0.2)];


IC_1 = [0.15;-.5];
IC_2 = [0.15;-.45];
[tall_1, Yall_1] = ode45(f,[0,100], IC_1);
[tall_2, Yall_2] = ode45(f,[0,100], IC_2);
figure
box on
plot(Yall_1(:,1), Yall_1(:,2), 'red', 'LineWidth', 1)
hold on
plot(Yall_2(:,1), Yall_2(:,2), 'blue', 'LineWidth', 1)
hold on
vectfieldn(f,-1:.15:2,-2:.15:2)
hold on
fplot(@(x) x.*(1-x).*(x-2),[-.5,2], 'k', 'LineWidth',1)
hold on
plot([0.2 0.2],[-2 2], 'k', 'LineWidth', 1)
xlabel('x');
ylabel('y');
legend('initial condition of (0.15,-0.5)','initial condition of (0.15,-0.45)')

hold off

figure
plot(tall_1,Yall_1(:,1),'red')
hold on
plot(tall_1,Yall_1(:,2),'red.')
hold on
plot(tall_2,Yall_2(:,1),'blue')
hold on
plot(tall_2,Yall_2(:,2),'blue.')
xlabel('t');
ylabel('x,y');
legend('x for initial condition of (0.15,-0.5)','y for initial condition of (0.15,-0.5)','x for initial condition of (0.15,-0.45)','y for initial condition of (0.15,-0.45)')

grid on
% plot([0.2 0.2],[-2 2], 'k', 'LineWidth', 1)