%make the color map for van der pol flow kick

vdp = vanderpol_repeatedkickscalculationJune2023;

%% vertical kicks

tvalues1 = 10:0.1:30;
kvalues1 = -0.22:0.005:-0.18;
xvalues_vertical = vdp.colors(tvalues1,kvalues1,2,400)

%% make vertical kick picture
load('vdp_tau_kick_xvalues.mat')
[T,K] = meshgrid(tvalues,kvalues);
%contourf(T,-K,xvalues',50,'LineColor','none')
%clim([1-1/sqrt(3) 1+1/sqrt(3)])
% or use 
imagesc(tvalues,kvalues,xvalues')
clim([1-1/sqrt(3) 1+1/sqrt(3)])


%% %horiztonal kicks in two pieces
% tau from 0.1 to 5, k = 0 to .75
% tau from 5 to 30, k = .5 to .75
tvalues1 = [0.1:0.1:5];
kvalues1 = [0:0.02:.75];
xvalues1 = vdp.colors(tvalues1,kvalues1,1,400);

tvalues2 = [5:0.2:30];
kvalues2 = [0.4:0.02:0.75];
xvalues2 = vdp.colors(tvalues2,kvalues2,1,400);

%% save if you like it!
save('vdp_tau_kick_xvalues_horiztonal_kick.mat','tvalues1','tvalues2','kvalues1','kvalues2','xvalues1','xvalues2')
%% actually plotting
%load('vdp_tau_kick_xvalues_horiztonal_kick.mat')
[T1,K1] = meshgrid(tvalues1,kvalues1);
[T2,K2] = meshgrid(tvalues2,kvalues2);
contourf(T1,K1,xvalues1',50,'LineColor','none')
clim([1-1/sqrt(3) 1+1/sqrt(3)])
hold on
contourf(T2,K2,xvalues2',50,'LineColor','none')
clim([1-1/sqrt(3) 1+1/sqrt(3)])