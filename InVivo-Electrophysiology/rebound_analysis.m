function [psth, delta] = rebound_analysis(file, bins, ind)
% Function to for each cell determine spike locations upon inactivation of
% a 590 nm light stimulus and to calculate the number of spikes in each bin 
% to create a peri-stimulus time histogram starting 100 ms before
% inactivation occurs and ending 100 ms after
%
% Inputs:
%     file - file that contains data recorded from the cell
%     lead - lead time prior to the stimulation used for segments
%     bins - number of bins for the peri-stimulus time histogram
%     artifact_template - template of the artifact for thalamic stimulation
%     ind - indices that indicate when thalamic stimulations occured
%
% Outputs:
%     psth - number of spikes per bin for the peri-stimulus time histogram
%     psth_spikeheight - average spikeheight for each bin
%     spikelocs - location of spikes to create a raster plot
%     si - sampling interval in us
%     delta - bin width in secs

% Initiate variables to store information for peri-stimulus time histogram
psth = zeros(1, bins);
psth_spikeheight = zeros(1, bins);

% Load the file and get the size of the data
[d,si,h] = abfload(file);
[p,q,r] = size(d);

% Determine the number of data points equivalent to 100 ms
len = 0.1/(si/1e6);

% Initiate variable to store data points from 100 ms before till 100 ms
% after inactivation of the 590 nm light stimulus
data = zeros(2*len, r);

% Extract and store the data
for i = 1:r
    data(:, r) = d(ind(2)-len+1:ind(2)+len, i);
end

% Get the size of the data
[p,r] = size(data);

% Calculate the size of each bin
delta = (p*(si/1e6))/bins;

% Initiate a variable to store the normalized data
data_normalized = zeros(p,r);

% Initiate a cell array to store the location of spikes for each
% inactivation of the 590 nm light stimulus
spikelocs = cell(1,q*r);

% Normalize data by subtracting the mean activity across the segment
for i = 1:r
    data_normalized(:,i) = data(:,i) - mean(data(:,i));
end

% For each segment find spike locations, store the in a cell array, 
% increase the number of spikes per bin for each spike that falls into 
% a specific bin and add the spikeheight for each spike to the correct bin
for i = 1:r
    [spikeheight, spikeloc] = findpeaks(data_normalized(:,i),'MinPeakProminence',0.6,'MaxPeakWidth',0.006/(si/1e6));
    spikelocs{i} = spikeloc;
    for j = 1:length(spikeloc)
        ind_bin = ceil(spikeloc(j)/(delta/(si/1e6)));
        psth(ind_bin) = psth(ind_bin) + 1;
        psth_spikeheight(ind_bin) = psth_spikeheight(ind_bin) + spikeheight(j);
    end
end

% Calculate the average spike height
psth_spikeheight = psth_spikeheight./psth;
psth_spikeheight(isnan(psth_spikeheight)) = 0;

% Normalize the number of spikes per bin by dividing through N*delta
psth = psth./(r*delta*1e3);

end