function[] = toadProc(directory)

%TOADPROC4 Takes directory of toad data files and saves tongue error data and
%figures

%% Find index of toad movies
if ~strcmp(pwd,directory)
    cd(directory);
end

if strcmp(computer,'MACI64')
    toadDir = cell2mat(strsplit(ls('-d','*b*cam2*'))');
    screenDir = deblank(ls('-d','*cal*_image*cam2*'));
%     toadDir = ls('-d',[pwd filesep '*b*cam2*']);
%     screenDir = deblank(ls('-d',[pwd filesep '*cal*_image*cam2*']));
else
    toadDir = ls([pwd filesep '*b*cam2*']);
    screenDir = ls([pwd filesep '*cal*_image*cam2*']);
end

    
%% run analysis

load([directory filesep screenDir filesep 'plane.mat']);
% screenData = importdata([directory filesep screenDir filesep 'DLTdv5_data_xyzpts.csv']);
for i = 1:size(toadDir,1)
    toad_folder = deblank(toadDir(i,:));
    fprintf('%s\n',['Processing ' toad_folder]);
    if ~exist([directory filesep toad_folder filesep 'DLTdv5_data_xyzpts.csv'])
        fprintf('%s\n','No xyz points file found');
        continue
    end
    toadData = importdata([directory filesep toad_folder filesep 'DLTdv5_data_xyzpts.csv']);
    if strcmp(computer,'MACI64')
        load([directory filesep deblank(ls([toad_folder filesep 'stim parameters*']))]); 
    else
        load([directory filesep toad_folder filesep ls([toad_folder filesep 'stim parameters*'])]);
    end
    [metrics] = tongueError(toadData,plane,touchData,params,toad_folder);
    save([toad_folder filesep 'metrics.mat'],'metrics');
    close all
    fprintf('%s\n','Processing complete!');
end