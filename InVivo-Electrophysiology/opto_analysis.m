function [psth, psth_spikeheight, spikelocs, si, delta] = opto_analysis(file, lead, bins, artifact_template, ind)
% Function to for each cell determine spike locations for each thalamic
% stimulation delivered, when thalamic stimulations are delivered in
% conjunction with a 590 nm light stimulus and to calculate the number
% of spikes in each bin to create a peri-stimulus time histogram
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

% Get the length from the beginning of one stimulation to the
% beginning of the next using the provided indices
trial_len = ind(3)-ind(1);

% Initiate variables necessary to reshape data such that
% there is data of one stimulation along the X-axis and stimulations 
% from one set of 10 stimulations along the Y axis
data = zeros(trial_len, length(ind)/10, r*5);
ind = vertcat(ind,p);
m = 1;

% Reshape data such that there is data of one stimulation along the X-axis
% and stimulations from one set of 10 stimulations along the Y axis
for i = 1:r
    k = 1;
    for j = 1:2:length(ind)-1
        if (ind(j+2)-ind(j)) == trial_len
           data(:, k, m) = d(ind(j)-lead/(si/1e6)+1:ind(j+2)-lead/(si/1e6), 1, i);
           k = k+1;
        else
           data(:, k, m) = d(ind(j)-lead/(si/1e6)+1:ind(j)-lead/(si/1e6)+trial_len, 1, i);
           k = 1;
           m = m+1;
        end
    end
end

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