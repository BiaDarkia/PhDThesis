function [psth, psth_spikeheight, spikelocs, si, delta] = mtstim_analysis(file, lead, bins, artifact_template, ind)
% Function to for each cell determine spike locations for each thalamic
% stimulation delivered and to calculate the number of spikes in each bin
% to create a peri-stimulus time histogram
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

% Reshape the data into segments such that there is data of one
% stimulation along the X-axis, stimulations from one sweep along
% the Y axis, and data from sweeps along the Z-axis
data = reshape(d(ind(1)-lead/(si/1e6)+1:ind(end-1)-lead/(si/1e6)+(ind(3)-ind(1)),1,:),[],10,r);

% Get the size of the data
[p,q,r] = size(data);

% Calculate the size of each bin
delta = (p*(si/1e6))/bins;

% Calculate starting point to use for normalizing data and
% initiate a variable to store the normalized data
normalize = p - 0.04/(si/1e6);
data_normalized = zeros(p,q,r);

% Initiate a cell array to store the location of spikes in response
% to each thalamic stimulation
spikelocs = cell(1,q*r);

% Normalize data by subtracting the mean activity across the last
% 40 ms of each segment as well as the artifact template
for i = 1:q
    for j = 1:r
        data_normalized(:,i,j) = data(:,i,j) - mean(data(normalize:end,i,j)) - artifact_template(:);
    end
end

% For each segment starting the lead prior to the thalamic ctimulation
% and lasting 100 ms - lead after the thalamic stimulation find
% spike locations, store the in a cell array, increase the number
% of spikes per bin for each spike that falls into a specific bin and
% add the spikeheight for each spike to the correct bin
for i = 1:r
    for j = 1:q
        ind_all = ((i-1)*q)+j;
        [spikeheight, spikeloc] = findpeaks(data_normalized(:,j,i),'MinPeakProminence',0.6,'MaxPeakWidth',0.006/(si/1e6));
        spikelocs{ind_all} = spikeloc;
        for k = 1:length(spikeloc)
            ind_bin = ceil(spikeloc(k)/(delta/(si/1e6)));
            psth(ind_bin) = psth(ind_bin) + 1;
            psth_spikeheight(ind_bin) = psth_spikeheight(ind_bin) + spikeheight(k);
        end
    end
end

% Calculate the average spike height
psth_spikeheight = psth_spikeheight./psth;
psth_spikeheight(isnan(psth_spikeheight)) = 0;

% Normalize the number of spikes per bin by dividing through N*delta
psth = psth./(q*r*delta*1e3);

end