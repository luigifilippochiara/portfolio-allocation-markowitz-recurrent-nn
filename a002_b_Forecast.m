%% 2.b Forecast
clear
clc
load data.mat

%% core

% bench
M=5000;     % number of simulated paths
N=60;      % length of simulated paths
w=6;        % average block length

% boostrapping row indices
[indices,~]=stationary_bootstrap_MC((1:size(returns,1))',M,w,N);

% compute ew portfolio
ewWeights = ones(1, size(returns,2))/size(returns,2);
ewReturns = returns*ewWeights';

benchSimReturns=zeros(N,M);
benchSimReturnsDefl=zeros(N,M);
ewSimReturns=zeros(N,M);
ewSimReturnsDefl=zeros(N,M);
for j=1:M
    inflationRate = inflation(indices(:,j));
    benchSimReturns(:,j)=cumprod(1+benchReturns(indices(:,j))./100);
    benchSimReturnsDefl(:,j)=cumprod(1+benchReturns(indices(:,j))./100-inflationRate)./(inflationRate+1);
    ewSimReturns(:,j)=cumprod(1+ewReturns(indices(:,j))./100);
    ewSimReturnsDefl(:,j)=cumprod(1+ewReturns(indices(:,j))./100-inflationRate)./(inflationRate+1);
end

period = dm.dates(end) + calmonths(1:N);
C = linspecer(2);

figure
% EW plot
subplot(1,2,1);
Q=quantile(ewSimReturns',[0.05 0.5 0.95])';
Qd=quantile(ewSimReturnsDefl',[0.05 0.5 0.95])';
fill([period, fliplr(period)], [Q(:,1)', fliplr(Q(:,3)')], C(1,:), 'FaceAlpha',.2, 'EdgeAlpha',0);
hold on
fill([period, fliplr(period)], [Qd(:,1)', fliplr(Qd(:,3)')], C(2,:), 'FaceAlpha',.2, 'EdgeAlpha',0);
plot(period,Q(:,2),'-','Color', C(1,:), 'LineWidth',2);
plot(period,Q(:,1),'--','Color', C(1,:),'LineWidth',2);
plot(period,Q(:,3),'--','Color', C(1,:),'LineWidth',2);
plot(period,Qd(:,2),'-','Color', C(2,:),'LineWidth',2);
plot(period,Qd(:,1),'--','Color', C(2,:),'LineWidth',2);
plot(period,Qd(:,3),'--','Color', C(2,:),'LineWidth',2);
title('5 years EW monthly returns');
ylabel('Return');
legend({'EW returns', 'Deflated EW returns'});
grid on

% EW plot
subplot(1,2,2);
Q=quantile(benchSimReturns',[0.05 0.5 0.95])';
Qd=quantile(benchSimReturnsDefl',[0.05 0.5 0.95])';
fill([period, fliplr(period)], [Q(:,1)', fliplr(Q(:,3)')], C(1,:), 'FaceAlpha',.2, 'EdgeAlpha',0);
hold on
fill([period, fliplr(period)], [Qd(:,1)', fliplr(Qd(:,3)')], C(2,:), 'FaceAlpha',.2, 'EdgeAlpha',0);
plot(period,Q(:,2),'-','Color', C(1,:), 'LineWidth',2);
plot(period,Q(:,1),'--','Color', C(1,:),'LineWidth',2);
plot(period,Q(:,3),'--','Color', C(1,:),'LineWidth',2);
plot(period,Qd(:,2),'-','Color', C(2,:),'LineWidth',2);
plot(period,Qd(:,1),'--','Color', C(2,:),'LineWidth',2);
plot(period,Qd(:,3),'--','Color', C(2,:),'LineWidth',2);
title('5 years Bench monthly returns');
ylabel('Return');
legend({'Bench returns', 'Deflated Bench returns'});
grid on

%% Probability of positive return
lastVal = ewSimReturns(end,:);
ewProb=length(lastVal(lastVal>1))/length(lastVal);
lastVal = ewSimReturnsDefl(end,:);
ewDeflProb=length(lastVal(lastVal>1))/length(lastVal);
lastVal = benchSimReturns(end,:);
benchProb=length(lastVal(lastVal>1))/length(lastVal);
lastVal = benchSimReturnsDefl(end,:);
benchDeflProb=length(lastVal(lastVal>1))/length(lastVal);

Tab1=table({'P(Returns > 1)'},ewProb,ewDeflProb,benchProb,benchDeflProb);
Tab1