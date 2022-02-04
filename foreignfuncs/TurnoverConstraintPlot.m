function TurnoverConstraintPlot(P, turnovers)
% TURNOVERCONSTRAINTPLOT: A helper function for generating several
% turnover-constrained portfolios on the same set of axes.  Inputs are P
% (the portfolio wihtout any turnover constaints) and turnovers (the
% desired thresholds)

% How smooth of efficient frontiers do we want?
numPorts = 10;

% An added benefit of the Portfolio object is the ability to easily change
% solvers.  Certain solvers are better suited for solving certain problems.
% In this case, switching from the default 'lcprog' solver to the
% interior-point-convex solver 'quadprog' gives us about a 450% speedup.
P = P.setSolver('quadprog');

% Generate frontier and moments for initial portfolio:
[Risk0, Return0]   = P.estimatePortMoments(P.InitPort);
scatter(Risk0, Return0, 40, 'g', 'filled', ...
    'DisplayName', 'Initial Portfolio')

hold all
% Generate unconstrained portfolio:
Wts = P.estimateFrontier(numPorts);
[Risks, Returns]   = P.estimatePortMoments(Wts);
plot(Risks, Returns, 'b', 'DisplayName', 'Efficient Frontier', ...
    'LineWidth', 2)

% Generate turnover-constrained portfolios
for currentTurnover = turnovers;
    P = P.setTurnover(currentTurnover);
    Wts = P.estimateFrontier(numPorts);
    [Risks, Returns]   = P.estimatePortMoments(Wts);
    name = [num2str(100*currentTurnover) '% Turnover'];
    plot(Risks, Returns, 'DisplayName', name, 'LineWidth', 2)
end
hold off
grid on
legend('Location', 'East')
title(P.Name)
xlabel('Standard Deviation of Portfolio Returns')
ylabel('Mean of Portfolio Returns')