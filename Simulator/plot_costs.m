function fig = plot_costs(J_step)

fig = figure;
set(gcf, 'color', 'w')
hold on;

yyaxis left
plot(cumsum(J_step), 'linewidth', 1);

axes = gca;
axes.XGrid = 'on';
axes.YGrid = 'off';
    
yyaxis right
plot(J_step, 'linewidth', 1);

xlim([1, length(J_step)]);
my_legend = legend('cumsum($J$, $t$)', '$J(t)$');
set(my_legend, 'Interpreter', 'latex');
xlabel('$t$', 'Interpreter', 'latex');