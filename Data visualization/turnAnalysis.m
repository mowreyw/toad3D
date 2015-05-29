function[] = turnAnalysis(directory)

%constants
traj_frames = 1250;

%load data
if ~strcmp(pwd,directory)
    cd(directory);
end

if strcmp(computer,'MACI64')
    toadDir = ls('-d',[pwd filesep '*b*cam2*']);
else
    toadDir = ls([pwd filesep '*b*cam2*']);
end

%start figure
n_recs = size(toadDir,1);
projection_time = NaN(n_recs,1);
tongue_error_uv = NaN(n_recs,2);
tongue_error_xyz = NaN(n_recs,3);
for i = 1:n_recs
    toad_folder = deblank(toadDir(i,:));
    if exist([directory filesep toad_folder filesep 'metrics.mat'],'file')
        load([directory filesep toad_folder filesep 'metrics.mat']);
    else
        continue
    end
    d_speed = diff([metrics.frame_vars.x metrics.frame_vars.y]);
    direction = sign(sum(d_speed(:,1)));
    speed = mean(sqrt(sum(d_speed.^2,2)));
%     speed = mean(hypot(d_speed(:,1),d_speed(:,2)));
    turn_time = find(diff(metrics.frame_vars.x));
    movie_turn_time = turn_time(end) - metrics.start_frame; %turn frame relative to start of recording
    projection_time(i) = metrics.t_hit - movie_turn_time;
    turn_position_uv = [metrics.frame_vars.x(turn_time(end))...
        metrics.frame_vars.y(turn_time(end))];
    turn_position_xyz = feval(metrics.plane.getXYZ,turn_position_uv(1),turn_position_uv(2));
%     turn_position_uv = metrics.target_coors_uv.center(turn_time,:);
%     turn_position_xyz = metrics.target_coors.center(turn_time,:);
    tongue_error_uv(i,:) = metrics.tongue_coors.uv(1:2)' - turn_position_uv;
    tongue_error_uv(i,1) = tongue_error_uv(i,1)*direction; %flip errors depending on heading direciton of target
    tongue_error_xyz(i,:) = metrics.tongue_coors.xyz - turn_position_xyz';
end

f1 = figure;
traj = (-speed*traj_frames:speed:0);
traj_uv = [ [traj zeros(1,traj_frames)] ; [zeros(1,traj_frames) -fliplr(traj)] ];
scatter(traj_uv(1,:),traj_uv(2,:),'CData',(-traj_frames:1:traj_frames));
hold on
scatter(tongue_error_uv(:,1),tongue_error_uv(:,2),'x','CData',projection_time);
caxis([-traj_frames traj_frames]);
colormap cool
axis equal
xlabel('u')
ylabel('v')


