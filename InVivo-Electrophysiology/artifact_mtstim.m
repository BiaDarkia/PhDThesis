function [artifact_template] = artifact_mtstim(file, lead, ind)
% Function to extract a template of the artifact produced by thalamic
% stimulation with a specific stimulation intensity and duration
% from a file that recorded this artifact in the absence of any cells.
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

% Reshape the data into segments such that there is data of one
% stimulation along the X-axis, stimulations from one sweep along
% the Y axis, and data from sweeps along the Z-axis
data = reshape(d(ind(1)-lead/(si/1e6)+1:ind(end-1)-lead/(si/1e6)+(ind(3)-ind(1)),1,:),[],10,r);   

% Average the artifact across all hundred stimulations in the file
% to create a template of the artifact
artifact_template = mean(data, [2 3]);

end