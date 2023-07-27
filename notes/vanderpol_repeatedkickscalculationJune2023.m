%% Van der pol 
% Use an initial virus level set to the kick size
% Kick repeatedly after time tau
% Look for an excursion (or not) at any point within numkick iterates

function funcs = vanderpol_repeatedkickscalculationJune2023
% Call this function first
% This gives access to all the functions in this file from the command line.
% Call with syntax of funcs = vanderpol_repeatedkickscalculationJune2023;
funcs.flow = @flow;
funcs.kickflowexcursion = @kickflowexcursion;
funcs.plot = @plotflowkick;
funcs.colors = @threshold_colors;
funcs.stkicks = @stochastic_kicks;
end

function [outputs] = stochastic_kicks(taurange,krange,kickd)
% didn't finish writing this :(

end

function [xvec] = threshold_colors(timevector,kickinputs,kickd,maxtime)
% kickd is kick direction. 1 for horiztonal, 2 for vertical
% 
xvec=[];
for k1=1:length(timevector)
    for j1=1:length(kickinputs)
        tau = timevector(k1);
        kick = kickinputs(j1);
        numkick = ceil(maxtime/tau);
        xvec(k1,j1) = kickflowexcursion(kick, kickd, tau, numkick);
        
    end
    timevector(k1)
end
end

function plotflowkick(kickv,kickd,tau,numkick)
% Plot flows and kicks.  Use plots to check the calculations.
% kickd should be 1 for horizontal and 2 for vertical
allt = [0];
allY = [];

%kick from equilibrium point
IC = [0.2; -0.288];

% Set initial virus to kick size:
IC(kickd) = IC(kickd)+kickv;
%tall = [tall; ts+sum(tau(1:(c-1)))];

% Kick repeatedly
for k1 = 1:numkick
    [t1,Y1] = flow(tau,IC);
    allt = [allt;t1+allt(end)];
    allY = [allY; Y1];
    IC = Y1(end,:);
    IC(kickd) = IC(kickd) + kickv;
end
allt = allt(2:end);
% Plot
figure
plot(allt, allY(:,1))
hold on 
plot(allt,allY(:,2))
xlabel('time')
legend('x','y')
title(['k = ', num2str(kickv), 'tau = ',num2str(tau)])
hold off
figure
plot(allY(:,1),allY(:,2))
hold
xlabel('x')
ylabel('y')
fplot(@(x) x.*(1-x).*(x-2),[-1,3],'k')
xline(0.2,'k')
title(['k = ', num2str(kickv), 'tau = ',num2str(tau)])
hold off
end


function [max_x] = kickflowexcursion(kick, kickd, tau, numkick)
%kickd is kick direction. 1 for horizontal or 2 for vertical.
% Initial virus and initial excursion:
allt = [0];
allY = [];

%kick from equilibrium point
IC = [0.2; -0.288];

% Set initial virus to kick size:
IC(kickd) = IC(kickd)+kick;
%tall = [tall; ts+sum(tau(1:(c-1)))];

% Kick repeatedly
for k1 = 1:numkick
    [t1,Y1] = flow(tau,IC);
    allt = [allt;t1+allt(end)];
    allY = [allY; Y1];
    IC = Y1(end,:);
    IC(kickd) = IC(kickd) + kick;
end
allt = allt(2:end);
max_x = max(allY(:,1));
end


function dydtval = dydt(t,Y)
%van der pol system

% variable order is 
%  [virus,target cells,infected cells,interferon,B cells,antibodies]
dydtval = [Y(1)*(1-Y(1))*(Y(1)-2)-Y(2);...
    0.1*(Y(1)-0.2)];
end

function [ts,YS] = flow(tau,IC)
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
[ts,YS] = ode15s(@(t,y) dydt(t,y),[0,tau],IC,options);
end