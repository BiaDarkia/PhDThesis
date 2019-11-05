function mtstim_plot(file, sweep, stim, ylim1, ylim2, lead, artifact_template, ind)
% Function to plot a data from a single sweep and a single thalamic stimulation.
%
% Inputs:
%     file - file that contains the data to be plotted
%     sweep - sweep of the data to be plotted
%     stim - number of thalamic stimulation in the sweep of the data to be plotted
%     ylim1 - lower limit of the y-axis for plotting
%     ylim 2 - upper limit of the y-axis for plotting
%     lead - lead time prior to the stimulation used for segments
%     artifact_template - template of the artifact for thalamic stimulation
%     ind - indices that indicate when thalamic stimulations occured

% Load the file and get the size of the data
[d,si,h] = abfload(file);
[p,q,r] = size(d);

% Plot data of a 2-second-long segment of activity at the end of the sweep between delivery of thalamic stimulations
data_long = d(2/(si/1e6):4/(si/1e6)-1, 1, sweep);

fig_name = strrep(strrep(file,'mtstim/',''),'.abf','');

data_long_fig = figure('visible','off','PaperPosition',[0 0 8 4]);
plot(1:500000, data_long)
hold on
xlim([0, 2/(si/1e6)])
xlabel('Time (s)')
xticks([0 0.5 1 1.5 2]./(si/1e6))
xticklabels({'0', '0.5', '1', '1.5', '2'})
ylabel('Activity (mV)')
hold off

% Save the plot
print(data_long_fig, strcat('figures/figure_', fig_name, '_long_', int2str(sweep)), '-dtiff', '-r300')


% Reshape the data into segments such that there is data of one
% stimulation along the X-axis, stimulations from one sweep along
% the Y axis, and data from sweeps along the Z-axis
data = reshape(d(ind(1)-lead/(si/1e6)+1:ind(end-1)-lead/(si/1e6)+(ind(3)-ind(1)),1,:),[],10,r);

% Get the size of the data
[p,q,r] = size(data);

% Calculate starting point to use for normalizing data and
% initiate a variable to store the normalized data
normalize = p - 0.04/(si/1e6);
data_normalized = zeros(p,q,r);

% Normalize data by subtracting the mean activity across the last
% 40 ms of each segment as well as the artifact template
for i = 1:q
    for j = 1:r
        data_normalized(:,i,j) = data(:,i,j) - mean(data(normalize:end,i,j)) - artifact_template;
    end
end

% Use the findepeaks algorithm to find spike locations
[spikeheight, spikeloc] = findpeaks(data_normalized(:,stim,sweep),'MinPeakProminence',0.6,'MaxPeakWidth',0.006/(si/1e6));

% Determine locations for x-ticks and figure name
tick_locs = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]./(si/1e6);

% Plot the data and indicate each spike location with a short black line
mtstim_fig = figure('visible','off','PaperPosition',[0 0 5 4]);
plot((1:p), data_normalized(:,stim,sweep), 'color', 'b')
hold on
xlim([0.01/(si/1e6), 0.09/(si/1e6)])
ylim([ylim1, ylim2])
for loc = spikeloc.'
    line([loc loc],[ylim2-1.5 ylim2-0.5], 'color', 'k', 'LineWidth', 2)
end
rectangle('Position',[0.019/(si/1e6) ylim1 0.002/(si/1e6) ylim2-ylim1], 'FaceColor', 'white', 'LineStyle', 'none')
xlabel('Time (ms)')
xticks(tick_locs)
xticklabels({'-10', '0', '10', '20', '30', '40', '50', '60', '70'})
ylabel('Activity (mV)')
hold off

% Save the plot
print(mtstim_fig, strcat('figures/figure_', fig_name, '_', int2str(sweep), '_', int2str(stim)), '-dtiff', '-r300')

