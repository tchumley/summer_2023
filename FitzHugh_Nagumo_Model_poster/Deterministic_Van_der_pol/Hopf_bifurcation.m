f=@(t,Y) [-Y(1)*(Y(1)-1)*(Y(1)-2)-Y(2);...
    0.1*(Y(1)-0.2)-0.03];
% (1\sqrt(3))



IC_1 = [1.13;0.02];
IC_2 = [2;-0.3];
[tall_1, Yall_1] = ode45(f,[0,900], IC_1);
[tall_2, Yall_2] = ode45(f,[0,1000], IC_2);
figure
box on
plot(Yall_1(:,1), Yall_1(:,2), 'red', 'LineWidth', 1)
hold on
plot(Yall_2(:,1), Yall_2(:,2), 'blue', 'LineWidth', 1)
hold on
vectfieldn(f,-1:.15:2.5,-2:.15:2)
hold on
fplot(@(x) x.*(1-x).*(x-2), 'k', 'LineWidth',1)
axis([-1 2.5 -2 2])
% hold on
plot([1-1/sqrt(3) 1-1/sqrt(3)],[-2 2], 'k', 'LineWidth', 1)
xlabel('x');
ylabel('y');
legend('initial condition of (1.13,0.02)','initial condition of (2,-0.3)')

hold off

% figure
% plot(tall_1,Yall_1(:,1),'red')
% hold on
% plot(tall_1,Yall_1(:,2),'red.')
% hold on
% plot(tall_2,Yall_2(:,1),'blue')
% hold on
% plot(tall_2,Yall_2(:,2),'blue.')
% xlabel('t');
% ylabel('x,y');
% legend('x for initial condition of (0.15,-0.5)','y for initial condition of (0.15,-0.5)','x for initial condition of (0.15,-0.45)','y for initial condition of (0.15,-0.45)')
% 
% grid on
% plot([0.2 0.2],[-2 2], 'k', 'LineWidth', 1)