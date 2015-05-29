function[] = reversalAnalysis(directory)

%constants
traj_frames = 800;

%load data
if ~strcmp(pwd,directory)
    cd(directory);
end
d = dir;

f1 = figure;
f2 = figure;
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
    
    %start figure
    n_recs = size(toadDir,1);
    projection_time = NaN(n_recs,1);
    tongue_error_uv = NaN(n_recs,1);
    tongue_proj = tongue_error_uv;
    % tongue_error_xyz = NaN(n_recs,3);
    toad_id = [];
    ang_speed = NaN(n_recs,1);
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
            toad_id{end + 1} = toad_folder;
            fprintf('Processing %s\n',toad_folder);
        else
            continue
        end
        d_speed = diff(metrics.frame_vars.x);
        d_speed = [d_speed(1); d_speed];
%         direction = sign(d_speed(1,1));
        direction = sign(d_speed);
        turn_time = find(diff(d_speed));
        if isempty(turn_time)
           turn_time = length(metrics.frame_vars.x);
        else
           turn_time = turn_time(1);
        end
        
        speed = abs(mean(d_speed(1:turn_time)));
        %     speed = mean(sqrt(sum(d_speed.^2,2)));
        %     speed = mean(hypot(d_speed(:,1),d_speed(:,2)));
        %     turn_time = find(diff(metrics.frame_vars.x));
        %     movie_turn_time = turn_time(end) - metrics.start_frame; %turn frame relative to start of recording
%         movie_turn_time = turn_time - metrics.start_frame; %start frane does not work correctly
%         d_target = diff(metrics.target_coors_uv.center(:,1));
%         target_moves = find(d_target);
%         target_start = target_moves(1);
        projection_time(i) = (metrics.t_hit - (turn_time +...
            metrics.ftrackbox.first_valid_frame - metrics.start_frame)); 
%         projection_time(i) = metrics.t_hit - (turn_time +...
%             metrics.ftrackbox.first_valid_frame); 
%         projection_time(i) = metrics.t_hit - movie_turn_time;
%         turn_position_uv = metrics.frame_vars.x(turn_time);
        %     turn_position_uv = [metrics.frame_vars.x(turn_time(end))...
        %         metrics.frame_vars.y(turn_time(end))];
        %     turn_position_xyz = feval(metrics.plane.getXYZ,turn_position_uv(1),turn_position_uv(2));
        %     turn_position_uv = metrics.target_coors_uv.center(turn_time,:);
        %     turn_position_xyz = metrics.target_coors.center(turn_time,:);
%         tongue_error_uv(i) = metrics.tongue_coors.uv(1) - turn_position_uv;
        tongue_error_uv(i) = (metrics.tongue_coors.uv(1) -...
            metrics.target_coors_uv.center(metrics.t_hit)) * direction(2);
        tongue_proj(i) = tongue_error_uv(i) + abs(projection_time(i))*-speed;
%         tongue_error_uv(i) = tongue_error_uv(i,1)*direction(2); %flip errors depending on heading direciton of target
        %     tongue_error_xyz(i,:) = metrics.tongue_coors.xyz - turn_position_xyz';
%         figure(f2);
%         hold off
%         plot(metrics.frame_vars.x(metrics.start_frame:end),'b');
%         hold on
% %         plot(turn_time,metrics.frame_vars.x(turn_time),'ro')
%         plot(metrics.target_coors_uv.center(...
%             metrics.ftrackbox.first_valid_frame:end,1),'r')
%         plot(metrics.t_hit-metrics.ftrackbox.first_valid_frame,...
%             metrics.tongue_coors.uv(1),'rx');
%         ylim([0 1000])
%         axis equal
%         waitforbuttonpress;
    end
    
    
    if ~exist('metrics','var')
        continue
    end
    % traj_uv = [ [traj zeros(1,traj_frames)] ; [zeros(1,traj_frames) -fliplr(traj)] ];
    % scatter(traj_uv(1,:),traj_uv(2,:),'CData',(-traj_frames:1:traj_frames));
    
    
    % scatter(tongue_error_uv(:,1),tongue_error_uv(:,2),'x','CData',projection_time);
    % caxis([-traj_frames traj_frames]);
    % colormap cool
    figure(f1);
    scatter(projection_time./metrics.fps,tongue_proj,44,'.');
    hold on
    
end
figure(f1);
plot([-traj_frames 0 traj_frames]/metrics.fps,speed*[-traj_frames 0 -traj_frames],'k');
% traj = [(-speed*traj_frames:speed:0),(-speed:-speed:-speed*traj_frames)];
% scatter((-traj_frames:1:traj_frames),traj,'k.');
xlabel('Time (s)')
ylabel('X (px)')
