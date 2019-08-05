% Script to analyse in-vivo electrophysiology data obtained in anesthetized
% animals that compares the amount of spikes induced in prelimbic pyramidal
% neurons upon thalamic stimulation without (files starting with
% mtstim-...) and with (files starting with opto-...) delivery of a 590 nm
% light stimulus to prelimbic cortical layer 1.
% In addition, this script will create a peri-stimulus time histogram for
% XXX ms before and after the 590 nm light stimulus is inactivated to
% confirm that no rebound spiking occurred.

% Load indices of the times delivery of thalamic stimulation begins and
% ends. These times were previously extracted from one mtstim-... and one
% opto-... file and stored in a csv file, since timing of thalamic 
% stimulations does not change across recorded cells or animals
mtstim_ind = csvread('mtstim-neg50mA-05ms-stimulus-waveform-indices.csv');
opto_ind = csvread('opto-neg50mA-05ms-stimulus-waveform-indices.csv');

% Load indices of the times delivery of the 590 nm light stimulus begins and
% ends. These times were previously extracted from one mtstim-... and one
% opto-... file and stored in a csv file, since timing of the 590 nm light
% stimulus does not change across recorded cells or animals
rebound_ind = csvread('opto-neg50mA-05ms-opto-waveform-indices.csv');

% Open any mtstim-... files to be included in the data analysis
mtstim_files = dir('mtstim/*.abf');
mtstim_files_names = {mtstim_files(:).name};

% Open the matching opto-... files
opto_files = dir('opto/*.abf');
opto_files_names = {opto_files(:).name};

% Open files that contain a recording of the artefact produced by
% thalamic stimulation in the absence of any cells to create a template for
% the artefact
artifact_files = dir('artifact/*.abf');
artifact_files_names = {artifact_files(:).name};

% Create a variable to store artefact templates
artifact_database = {'index', 'mtstim_template', 'opto_template'};

% Specify the lead time prior to thalamic stimulation that will be included
% in each segment. Each segement will run from the -lead before until 
% 100 ms - lead after the thalamic stimulation.
lead = 0.02;

% Specify the number of bins for the peri-stimulus time histogram to
% analyse spiking activity after the thalamic stimulation
bins = 1000;

% Specify the number of bins for the peri-stimulus time histogram to
% analyse for rebound spiking upon inactivation of the 590 nm light
% stimulus
rebound_bins = 200;

% Data will be truncated around the thalamic stimulation. The trauncation
% starts at the beginning of the thalamic stimulation, which is equal to lead 
% from the start of each segment, and ends at trunc that is lead + XX secs
trunc = lead + 0.002;

% Data from each cell will be analyzed automatically using the functions
% mtstim_analysis, opto_analysis and rebound_analysis. Each of these 
% return data to create a normalized peri-stimulus time histogram. Here,
% variables to store and combine this data are initialized
mtstim_psth_total = zeros(1, bins);
opto_psth_total = zeros(1, bins);
rebound_psth_total = zeros(1, rebound_bins);

% Functions also return data on the spikeheight of spikes in each bin of
% the peri-stimulus time histogram. Here, variables to store and combine
% those are initiated
mtstim_psth_spikeheight_total = zeros(1, bins);
opto_psth_spikeheight_total = zeros(1, bins);

% Initiate necessary index and counting variables
ind_mtstim = 0;
ind_opto = 0;
mtstim_n = 0;
opto_n = 0;
n = 0;

% Initiate variables that will store spike locations for each cell
% to later create a raster plot
mtstim_raster = 0;
opto_raster = 0;

% For each cell open the mtstim-... file
for f = 1:length(mtstim_files_names)
    mtstim_file = mtstim_files_names{f};
    
    % Check if there already exists an artefact template in the database
    ind_artifact = find(contains(artifact_database(:, 1), string(mtstim_file(1:5))));
    if ~isempty(ind_artifact)
        artifact_template_mtstim = artifact_database{ind_artifact, 2};
        artifact_template_opto = artifact_database{ind_artifact, 3};
        
    % If there is none call the functions artifact_mtstim and artifact_opto
    % to extract an artefact template from the files containing the
    % recording of the artefact after thalamic stimulation in the presence
    % and absence of a 590 nm light stimulus, and save them to the artefact
    % template database
    else     
        ind_artifact_mtstim_pattern = strcat(mtstim_file(1:5), '[0-9_]+', 'artifact-mtstim');
        ind_artifact_mtstim = find(~cellfun(@isempty,regexp(artifact_files_names, ind_artifact_mtstim_pattern)));        
        ind_artifact_opto_pattern = strcat(mtstim_file(1:5), '[0-9_]+', 'artifact-opto');
        ind_artifact_opto = find(~cellfun(@isempty,regexp(artifact_files_names, ind_artifact_opto_pattern)));
        [artifact_template_mtstim] = artifact_mtstim(strcat('artifact/', artifact_files_names{ind_artifact_mtstim}), lead, mtstim_ind);        
        [artifact_template_opto] = artifact_opto(strcat('artifact/', artifact_files_names{ind_artifact_opto}), lead, opto_ind);
        [ind1, ind2] = size(artifact_database);
        artifact_database(ind1+1, :) = {mtstim_file(1:5), artifact_template_mtstim, artifact_template_opto};
    end
    
    % Find the opto-... file that matches the mtstim-... file, i.e.
    % contains recording data from the same cell
    ind_opto_pattern = strcat(mtstim_file(1:5), '[0-9_]+', mtstim_file(10:14));
    ind_opto_file = find(~cellfun(@isempty,regexp(opto_files_names, ind_opto_pattern)));
    
    % If there is a matching opto-... file, call the functions
    % mtstim_analysis, opto_analysis and rebound_analysis to extract spike
    % locations from the data and sort spikes into bins to later on create
    % peri-stimulus time histrograms and raster plots
    if ~isempty(ind_opto_file)
        opto_file = opto_files_names{ind_opto_file};
        
        [mtstim_psth, mtstim_psth_spikeheight, mtstim_spikelocs, mtstim_si, mtstim_delta] = mtstim_analysis(strcat('mtstim/', mtstim_file), lead, bins, artifact_template_mtstim, mtstim_ind);
        [opto_psth, opto_psth_spikeheight, opto_spikelocs, opto_si, opto_delta] = opto_analysis(strcat('opto/', opto_file), lead, bins, artifact_template_opto, opto_ind); 
        [rebound_psth, rebound_delta] = rebound_analysis(strcat('opto/', opto_file), rebound_bins, rebound_ind);
        
        % Generate a figure of a peri-stimulus time histogram that
        % summarizes spiking activity after thalamic stimulaion without and
        % with delivery of the 590 nm light stimulus
        fig_name = strcat('figures/psth-', mtstim_file(1:14));
        
        psth_fig = figure('visible','off');
        bar((-99:900)*mtstim_delta*1e3, mtstim_psth(1:1000))
        hold on
        bar((-99:900)*opto_delta*1e3, opto_psth(1:1000))
        ylim([0 1])
        hold off
        print(psth_fig, fig_name, '-dtiff', '-r300')
        
        % Add up data from each cell cluster into one array
        mtstim_psth_total = mtstim_psth_total + mtstim_psth;
        opto_psth_total = opto_psth_total + opto_psth;
        rebound_psth_total = rebound_psth_total + rebound_psth;
        mtstim_psth_spikeheight_total = mtstim_psth_spikeheight_total + mtstim_psth_spikeheight;
        opto_psth_spikeheight_total = opto_psth_spikeheight_total + opto_psth_spikeheight;
        
        % Store spike locations from each thalamic stimulation for that
        % cell into a matrix so that the first row holds the spike
        % locations and the second row holds a number as identifier for
        % each thalamic stiulation
        for i = 1:length(mtstim_spikelocs)
            if isempty(mtstim_spikelocs{i})
                ind_mtstim = ind_mtstim + 1;
                mtstim_raster(1, ind_mtstim) = -100;
                mtstim_raster(2, ind_mtstim) = mtstim_n + i;
            else
                for j = 1:length(mtstim_spikelocs{i})
                    ind_mtstim = ind_mtstim + 1;
                    mtstim_raster(1, ind_mtstim) = mtstim_spikelocs{i}(j)*(mtstim_si/(mtstim_delta*1e6));
                    mtstim_raster(2, ind_mtstim) = mtstim_n + i;
                end
            end
        end
        
        % Do the same for data extracted from recordings obtained upon
        % delivery of the 590 nm light stimulus
        for i = 1:length(opto_spikelocs)
            if isempty(opto_spikelocs{i})
                ind_opto = ind_opto + 1;
                opto_raster(1, ind_opto) = -100;
                opto_raster(2, ind_opto) = opto_n + i;
            else
                for j = 1:length(opto_spikelocs{i})
                    ind_opto = ind_opto + 1;
                    opto_raster(1, ind_opto) = opto_spikelocs{i}(j)*(opto_si/(opto_delta*1e6));
                    opto_raster(2, ind_opto) = opto_n + i;
                end
            end
        end
        
        % Advance counters by the number of spikes recorded from the cell
        mtstim_n = mtstim_n + length(mtstim_spikelocs);
        opto_n = opto_n + length(opto_spikelocs);
        
        % Advance counter by one for each cell
        n = n + 1;
        
    end
end

% Normalize peri-stimulus time histograms by dividing by the 
% number of cells recorded
mtstim_psth_total = mtstim_psth_total./n;
opto_psth_total = opto_psth_total./n;

% Truncate the artefact
mtstim_psth_total(lead/mtstim_delta:trunc/mtstim_delta) = 0;
opto_psth_total(lead/opto_delta:trunc/opto_delta) = 0;

% Plot an overlay of the peri-stimulus time histogram summarizing spiking
% activity slightly prior and after thalamic stimulation and save the figure
psth_fig = figure('visible','off','PaperPosition', [0 0 10 4]);
bar((-99:300)*mtstim_delta*1e3, mtstim_psth_total(101:500), 'b')
hold on
bar((-99:300)*opto_delta*1e3, opto_psth_total(101:500), 'r')
title({'A- Post-Stimulus Time Histogram after VM Stimulation',''})
xlabel('Time (ms)')
ylabel('Average number of spikes/ms')
legend('Light OFF', 'Light ON')
ylim([0 0.6])
hold off       
print(psth_fig, 'figures/figure_psth', '-dtiff', '-r300')

% Normalize the peri-stimulus time histogram for rebound spiking
% by dividing by the number of cells recorded
rebound_psth_total = rebound_psth_total./n;

% Calculate the mean spiking activity prior and after inactivation 
% of the 590 nm light stimulus
mean_prior_rebound = mean(rebound_psth_total(1:length(rebound_psth)/2));
mean_post_rebound = mean(rebound_psth_total((length(rebound_psth)/2):length(rebound_psth)));

% Plot the peri-stimulus time histogram summarizing spiking activity 
% 100 ms before and after inactivation of the 590 nm stimulus light
% and save the figure
psth_fig = figure('visible','off','PaperPosition', [0 0 10 4]);
bar((1:length(rebound_psth))*rebound_delta*1e3, rebound_psth_total, 'k')
hold on
line([1 length(rebound_psth)/2],[mean_prior_rebound mean_prior_rebound], 'color', 'red', 'LineStyle', '--', 'LineWidth', 1)
line([length(rebound_psth)/2 length(rebound_psth)],[mean_post_rebound mean_post_rebound], 'color', 'red', 'LineStyle', '--', 'LineWidth', 1)
line([length(rebound_psth)/2 length(rebound_psth)/2],[0 1], 'color', 'red', 'LineWidth', 2)
title({'C- Peri-Stimulus Time Histogram to Illustrate Absence of Rebound Spiking',''})
xlabel('Time (ms)')
xticks([0 50 100 150 200])
xticklabels({'-100','-50','0','50','100'})
ylabel('Average number of spikes/ms')
ylim([0 0.2])
hold off
print(psth_fig, 'figures/figure_psth_rebound', '-dtiff', '-r300')

% Normalize peri-stimulus time histograms summarizing the average spike 
% height by dividing by the number of cells recorded
mtstim_psth_spikeheight_total = mtstim_psth_spikeheight_total./n;
opto_psth_spikeheight_total = opto_psth_spikeheight_total./n;

% Truncate the artefact
mtstim_psth_spikeheight_total(lead/mtstim_delta:trunc/mtstim_delta) = 0;
opto_psth_spikeheight_total(lead/opto_delta:trunc/opto_delta) = 0;

% Plot an overlay of the peri-stimulus time histogram summarizing spike
% height slightly prior and after thalamic stimulation and save the figure
psth_fig = figure('visible','off','PaperPosition', [0 0 5 3]);
bar((-99:300)*mtstim_delta*1e3, mtstim_psth_spikeheight_total(101:500), 'b')
hold on
bar((-99:300)*opto_delta*1e3, opto_psth_spikeheight_total(101:500), 'r')
title({'Average Strength of Spikes after VM Stimulation',''})
xlabel('Time (ms)')
ylabel('Average Strength of Spikes (mV)')
legend('Light OFF', 'Light ON')
ylim([0 4])
hold off       
print(psth_fig, 'figures/figure_psth_spikeheight', '-dtiff', '-r300')

% Create two raster plots on top of each other to compare spiking activity
% a little bit prior to thalamic stimulation and after thalamic stimulation
% without and with the 590 nm light stimulus being activated
raster_fig = figure('visible','off','PaperPosition', [0 0 10 4]);
subplot(2,1,1);
scatter1 = scatter(mtstim_raster(1,:), mtstim_raster(2,:),'b','.');
hold on
rectangle('Position',[199 -40 21 2140], 'FaceColor', 'white', 'LineStyle', 'none')
xlabel('Time (ms)');
ylabel('Trials')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-10','-5', '0', '5', '10', '15', '20', '25', '30'})
ylim([-100 2100])
title({'B- Raster Plot to Illustrate Activity after VM Stimulation',''})
subplot(2,1,2);
scatter2 = scatter(opto_raster(1,:), opto_raster(2,:),'r','.');
hold on
rectangle('Position',[199 -40 21 2140], 'FaceColor', 'white', 'LineStyle', 'none')
xlabel('Time (ms)')
ylabel('Trials')
xlim([100 500])
xticks([100 150 200 250 300 350 400 450 500])
xticklabels({'-10','-5', '0', '5', '10', '15', '20', '25', '30'})
ylim([-100 2100])
legend([scatter2, scatter1],{'Light ON', 'Light OFF'},'Location','southeast','FontSize',6)
print(raster_fig, 'figures/figure_raster', '-dtiff', '-r300')

% Plot the recorded data from one cell for one thalamic stimulation
ind_artifact = find(contains(artifact_database(:, 1), '19717'));
artifact_template = artifact_database{ind_artifact, 2};
mtstim_plot('mtstim/19717025_28030-mtstim-neg50mA-05ms.abf', 6, 3, -2, 6, lead, artifact_template, mtstim_ind);

% Plot the recorded data from another cell for one thalamic stimulation
ind_artifact = find(contains(artifact_database(:, 1), '19728'));
artifact_template = artifact_database{ind_artifact, 2};
mtstim_plot('mtstim/19728012_21060-mtstim-neg50mA-05ms.abf', 5, 1, -15, 35, lead, artifact_template, mtstim_ind);

% Apply a Wilcoxon signed rank test to determine if the number of spikes
% per ms from 2 ms after thalamic stimulation, i.e. after the end of the 
% artefact, and till 20 ms after the end of the thalamic stimulation
% differs significantly when a 590 nm light stimulus is delivered as
% opposed to it not being delivered
p_spikes = signrank(mtstim_psth_total(221:400), opto_psth_total(221:400));

% Display the output of the test as well as the average number of spikes 
% per ms from 2 ms after thalamic stimulation, i.e. after the end of the 
% artefact, and till 20 ms after the end of the thalamic stimulation for
% the 590 nm light stimulus being delivered and not being delivered
disp(p_spikes)
disp(mean(mtstim_psth_total(221:400)))
disp(mean(opto_psth_total(221:400)))

% Display the average number of spikes per ms from 20 ms till 80 ms
% after the end of the thalamic stimulation for the 590 nm light stimulus
% being delivered and not being delivered
disp(mean(mtstim_psth_total(401:1000)))
disp(mean(opto_psth_total(401:1000)))
