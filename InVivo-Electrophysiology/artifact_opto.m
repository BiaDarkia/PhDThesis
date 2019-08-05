function [artifact_template] = artifact_opto(file, lead, ind)
% Function to extract a template of the artifact produced by thalamic
% stimulation with a specific stimulation intensity and duration upon
% exposure of the brain to a 590 nm light stimulus from a file 
% that recorded this artifact in the absence of any cells.
%
% Inputs:
%     file - file that contains the recording of the artifact
%     lead - lead time prior to the stimulation used for segments
%     ind - indices that indicate when thalamic stimulations occured
%
% Output:
%     artifact_template - template of the artifact

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

% Average the artifact across all hundred stimulations in the file
% to create a template of the artifact
artifact_template = mean(data, [2 3]);

end
