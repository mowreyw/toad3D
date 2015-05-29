function[mean_data,metrics_matrix] = meanMetric(directory,metric_name,dim)

%MEANMETRIC Calculates the mean of a metric in toad data file.
% Takes directory of toad data files and the name of a metric. Collects
% data from toad metric files and takes the mean. If metric is
% time-varying, takes the mean of values from some pre-turn interval. This
% leaves a row vector is metric is also multidimensional. Row vectors are
% concatenated into a matrix, and population mean is taken.

%Constants
pre_turn = 36;

%Collate data
if ~strcmp(pwd,directory)
    cd(directory);
end
d = dir;

metrics_matrix = [];
for j = 3:numel(d)
    if ~d(j).isdir
        continue
    end
    
    if strcmp(computer,'MACI64')
        toadDir = ls('-d',[d(j).name filesep '*b*cam2*']);
        toadDir = cell2mat(strsplit(toadDir)');
    else
        toadDir = ls([d(j).name filesep '*b*cam2*']);
    end
    
    if isempty(toadDir)
        continue
    end
    
    n_recs = size(toadDir,1);
%     metrics_matrix = [metrics_matrix; NaN(n_recs,dim)];
     
    for i = 1:n_recs
        toad_folder = deblank(toadDir(i,:));
        if strcmp(computer,'MACI64')
            fname = [directory filesep toad_folder filesep 'metrics.mat'];
        else
            fname = [directory filesep d(j).name filesep toad_folder...
                filesep 'metrics.mat'];
        end
        
        if exist(fname,'file')
            load(fname);
%             toad_id{end + 1} = toad_folder;
            fprintf('Processing %s\n',toad_folder);
        else
            continue
        end
        
        metrics_matrix = [metrics_matrix; nanmean(metrics.(metric_name)...
            (metrics.t_0-pre_turn:metrics.t_0,:))];    
    end
    mean_data = nanmean(metrics_matrix);
end
