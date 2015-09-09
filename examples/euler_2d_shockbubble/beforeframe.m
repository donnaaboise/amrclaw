figure(2);

plot_level = input('Input level to plot : ');
if (isempty(plot_level))
    PlotData = ones(1,MaxLevels);
else
    PlotData = [zeros(1,plot_level-1) 1 zeros(1,MaxLevels)];
end

