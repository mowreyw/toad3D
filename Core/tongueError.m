function[metrics]=tongueError(toadData,plane,touchData,params,toad_folder)

%TONGUEERROR3D Calculate metrics relevant to tongue projection error in 3D
%Takes as its arguement xyz points de csv files from Ty Hedricks. Plane is
%a structure with fields normal, basis, and point. Point should represent
%origin (center point) of plane.

%% constants

screen_size = [1280 720];
movie_frames = 400;
metrics = struct('centroid_euclidian_error',NaN,...
    'leading_edge_euclidian_error',NaN,'centroid_angular_error',NaN,...
    'leading_edge_angular_error',NaN,'alpha',NaN,'beta',NaN,'gamma',NaN,...
    'prey_azimuth',NaN,'prey_elevation',NaN,'prey_angular_size',NaN,...
    'prey_angular_velocity',NaN,'prey_angular_speed',NaN,...
    'prey_metric_size',NaN,'prey_metric_speed',NaN,'plane',plane);

ind = strfind(toad_folder,filesep);
if isempty(ind)
    fig_name = toad_folder;
elseif ind(end) == numel(toad_folder)
    fig_name = toad_folder(ind(end-1)+1:end-1);
else
    fig_name = toad_folder(ind(end)+1:end);
end
fig_name = strrep(fig_name,'_',' ');

if strcmp(params.fps_mode,'triple')
    fps = 360;
else
    fps = 120;
end
metrics.fps = fps;

%% load xyz datapoints

headPts = toadData.data(:,1:9);
metrics.headPts = headPts;
tonguePts = [];
if size(toadData.data,2) >= 15
    tonguePts = toadData.data(:,10:15);
    thit = find(~isnan(tonguePts(:,1)));
    headPts = headPts(1:thit,:);
    
    %find the tongue equation
    l_b = tonguePts(thit,1:3);
    l_a = tonguePts(thit,4:6);
    t_line = @(t) l_a + t*(l_b-l_a);
else
    thit = length(headPts);
end
metrics.t_hit = thit;
time_vec = [0:1:thit-1]./fps;


%% plot the intersection

fig1 = figure;
corners = [0 0; screen_size(1) 0; screen_size; 0 screen_size(2)];
corners_xyz = NaN(4,3);
for i = 1:4
   corners_xyz(i,:) = feval(plane.getXYZ,corners(i,1),corners(i,2))'; 
end
fill3(corners_xyz(:,1),corners_xyz(:,2),corners_xyz(:,3),zeros(4,1),...
    'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none');

hold on
fill3(headPts(thit,1:3:7),headPts(thit,2:3:8),headPts(thit,3:3:9),'k');
plot3(headPts(1:thit,1),headPts(1:thit,2),headPts(1:thit,3),'b')
plot3(headPts(1:thit,4),headPts(1:thit,5),headPts(1:thit,6),'r')
plot3(headPts(1:thit,7),headPts(1:thit,8),headPts(1:thit,9),'g')
% plot3(p_0(1),p_0(2),p_0(3),'bo')
% plot3(p_1(1),p_1(2),p_1(3),'ro')
% plot3(p_2(1),p_2(2),p_2(3),'go')
if ~isempty(tonguePts)
    t = ((plane.origin_xyz-l_a)*plane.normal) / ((l_b-l_a)*plane.normal);
    intersection = feval(t_line,t); %intersection with projection surface
    plot3(intersection(1),intersection(2),intersection(3),'rx')
    line([l_b(1) intersection(1)],[l_b(2) intersection(2)],[l_b(3) intersection(3)]);
end
set(gca,'ZDir','reverse','YDir','reverse','XGrid','On','YGrid','On',...
    'ZGrid','On')
xlabel('X'); ylabel('Y'); zlabel('Z');
% set(gca,'ZDir','reverse')
% set(gca,'XDir','reverse')
axis equal

%% filter photodiode data, register projector and camera frame counts
thresh = 2;
rec_end = find(touchData.record < thresh);
rec_end = rec_end(1); %end of camera recording

%photodiode signal
forder = 50;
cut = 0.01;
flt = fir1(forder,cut);
pd = filtfilt(flt,1,touchData.photodiode(1:rec_end));
pd = pd - min(pd(forder/2:end-forder/2)); %filtered and zeroed pd signal
p_thresh = max(pd(forder/2:end-forder/2))/2;
bin_pd = logical(pd > p_thresh); %binarized pd signal
dpd = [0; diff(bin_pd)];
pd_on = find(dpd(1:rec_end) > 0.5);
pd_ch = find(abs(dpd(1:rec_end)) > 0.5);
d_pd_ch = diff(pd_ch);
fprintf('Size of first gap %d\n',d_pd_ch(1));
if (d_pd_ch(1)) > (mean(d_pd_ch) + std(d_pd_ch))
    pd_ch(1) = []; %delete first pd signal. deals with pathological condition where pd is still high at start of recording.
end

%frame drops
rle = [find(bin_pd(2:end) - bin_pd(1:end-1)); length(bin_pd)]; %run-length encoding
runs = diff([0; rle]);
vals = bin_pd(rle);
norm_runs = round(runs ./ mode(runs));
drops = find(norm_runs > 1.5); %frame drop indices
drops = drops(2:end-1); %ignore pre-trial and post trial periods
n_drops = sum(norm_runs(drops)-1); %total number of dropped frames

%sync signal
sync = touchData.sync(1:rec_end) - min(touchData.sync);
sync_thresh = max(sync)/2;
bin_sync = logical(sync > sync_thresh);
dsync = [0; diff(bin_sync)];
sync_on = find(dsync > 0.5);
jitter = 10; %should be proportional to DAQ sampling rate, not just a constant!
index = sum(abs(repmat(pd_ch,1,numel(sync_on)) - repmat(sync_on',numel(pd_ch),1)) < jitter); %find frames where coincide within some jitter...
stim_moves = index' + circshift(index',1) + circshift(index',2);
% NB: not sure correct approach is being applied here. if dropped frames
% are never shown, then should probably reconstruct stimulus based on vsync
% signal and photodiode start and end.

%registration
record_start = sync_on(end-(movie_frames-1));
unrecorded_frames = sync_on(sync_on(logical(stim_moves)) < record_start);
recorded_frames = stim_moves(find(sync_on == record_start):end);
recorded_frames = recorded_frames(1:thit);
first_valid_frame = 1;
if pd_ch(1) > record_start
    first_valid_frame = sum((sync_on >= record_start) &...
        (sync_on < pd_ch(1))) + 1;
    recorded_frames(1:first_valid_frame-1) = NaN;
end

%save measures
metrics.ftrackbox.unrecorded_frames = unrecorded_frames;
metrics.ftrackbox.recorded_frames = recorded_frames;
metrics.ftrackbox.first_valid_frame = first_valid_frame;
metrics.ftrackbox.frame_drop_indices = drops;

%plot results
fig2 = figure;
plot(touchData.timeStamps,(touchData.sync - min(touchData.sync))/5,'r','LineWidth',2)
hold on
plot(touchData.timeStamps,(touchData.record - min(touchData.record))/5,'g','LineWidth',2)
plot(touchData.timeStamps(1:rec_end),bin_pd,'LineWidth',2)
plot(touchData.timeStamps(record_start),1,'yx','MarkerSize',10,'LineWidth',2)
legend('Camera Exposure','Camera record','Photodiode','Record start')
title([fig_name ' PD sync plot'])
saveas(fig2,[toad_folder filesep 'PD synchronization.fig']);

%find pixel position of target center and calculate error
if isempty(tonguePts)
    return
end

%% Plot target trajectory

if ~isfield(params,'frame_vars')
%     target_init = params.objXinit + (length(unrecorded_frames)*params.objVelX);
    target_init = [params.objXinit + (length(unrecorded_frames)*params.objVelX),...
        params.objYinit + (length(unrecorded_frames)*params.objVelY)];
%     traj_x = target_init+[NaN(first_valid_frame-1,1); cumsum(recorded_frames(first_valid_frame:end)*params.objVelX)];
    center = [target_init(1)+[NaN(first_valid_frame-1,1); cumsum(recorded_frames(first_valid_frame:end)*params.objVelX)],...
        target_init(2)+[NaN(first_valid_frame-1,1); cumsum(recorded_frames(first_valid_frame:end)*params.objVelY)]];
%     direction = sign(traj_x(end) - params.objXinit);
%     traj_lx = traj_x + direction*params.objLenX/2;
%     traj_tx = traj_x - direction*params.objLenX/2;
%     traj_y = ones(size(traj_x))*params.objYinit;
%     direction = [ 0 0 ; diff(center)./ repmat(sqrt(sum(diff(center).^2,2)),1,2)];
%     sq_dir = direction.^2;
    a = repmat(params.objLenX/2,length(unrecorded_frames) + length(recorded_frames),1);
    b = repmat(params.objLenY/2,length(unrecorded_frames) + length(recorded_frames),1);
%     coef = a.*b ./ sqrt(a.^2.*sq_dir(:,2)+ b.^2.*sq_dir(:,1));
%     leading_edge = center + repmat(coef,1,2).*direction;
%     trailing_edge = center - repmat(coef,1,2).*direction;
else
    fvars_file = dir([toad_folder filesep 'frame*']);
    T = readtable([toad_folder filesep fvars_file.name],'Delimiter',' ');
    start_frame = length(unrecorded_frames)+1;
%     traj_x = [NaN(first_valid_frame-1,1); T.x(start_frame+cumsum(recorded_frames(first_valid_frame:end)))];
%     traj_y = [NaN(first_valid_frame-1,1); T.y(start_frame+cumsum(recorded_frames(first_valid_frame:end)))];
    center = [[NaN(first_valid_frame-1,1); T.x(start_frame+cumsum(recorded_frames(first_valid_frame:end)))],...
        [NaN(first_valid_frame-1,1); T.y(start_frame+cumsum(recorded_frames(first_valid_frame:end)))]];
    a = [NaN(first_valid_frame-1,1); T.r1(start_frame+cumsum(recorded_frames(first_valid_frame:end)))/2];
    b = [NaN(first_valid_frame-1,1); T.r2(start_frame+cumsum(recorded_frames(first_valid_frame:end)))/2];
    metrics.frame_vars = T;
    metrics.start_frame = start_frame;
%     direction = [ 0 0 ; diff(center)./ repmat(sqrt(sum(diff(center).^2,2)),1,2)];
%     sq_dir = direction.^2;
%     coef = a.*b ./ sqrt(a.^2.*sq_dir(:,2)+ b.^2.*sq_dir(:,1));
%     leading_edge = center + repmat(coef,1,2).*direction;
%     trailing_edge = center - repmat(coef,1,2).*direction;
%     direction = sign(diff([traj_x traj_y]));
%     traj_lx = traj_x + params.objLenX/2;
%     traj_tx = traj_x - params.objLenX/2;
end
direction = [ 0 0 ; diff(center)./ repmat(sqrt(sum(diff(center).^2,2)),1,2)];
sq_dir = direction.^2;
coef = a.*b ./ sqrt(a.^2.*sq_dir(:,2)+ b.^2.*sq_dir(:,1));
leading_edge = center + repmat(coef,1,2).*direction;
trailing_edge = center - repmat(coef,1,2).*direction;
trajectory_frames = length(leading_edge);
top_edge = NaN(trajectory_frames,2);
bottom_edge = top_edge;

fn = {'center','leading_edge','trailing_edge','top_edge','bottom_edge'};
fd = { center; leading_edge; trailing_edge; top_edge; bottom_edge};
% fd = {traj_x, traj_y; traj_lx, traj_y; traj_tx, traj_y; traj_x, ...
%     traj_y + params.objLenY/2; traj_x, traj_y  - params.objLenY/2};

for i = 1:numel(fn)
    raw = NaN(trajectory_frames,3);
    warped = raw;
    for j = 1:trajectory_frames
        raw(j,:) = plane.uv2xyz*(([fd{i}(j,1); fd{i}(j,2); 1] - [plane.origin_uv 0]').*plane.mmpp');
        warped(j,:) = feval(plane.getXYZ,fd{i}(j,1),fd{i}(j,2))';
%         raw(j,:) = plane.uv2xyz*(([fd{i,1}(j); fd{i,2}(j); 1] - [plane.origin_uv 0]').*plane.mmpp');
%         warped(j,:) = feval(plane.getXYZ,fd{i,1}(j),fd{i,2}(j))';
    end
    metrics.target_coors_uv.(fn{i}) = fd{i};
    metrics.target_coors_raw.(fn{i}) = raw;
    metrics.target_coors.(fn{i}) = warped;
end
metrics.tongue_coors.xyz = intersection;
metrics.tongue_coors.uv = feval(plane.getUV,intersection');
% metrics.target_coors.uv = [traj_x, traj_y];

%Plot results
figure(fig1);
plot3(metrics.target_coors.center(:,1),metrics.target_coors.center(:,2),metrics.target_coors.center(:,3),'yo','MarkerSize',2)
plot3(metrics.target_coors.center(end,1),metrics.target_coors.center(end,2),metrics.target_coors.center(end,3),'kx','MarkerSize',12)
% plot3(metrics.touch_coors_3D(1),metrics.touch_coors_3D(2),metrics.touch_coors_3D(3),'x','MarkerSize',12,'MarkerEdgeColor',[0 0.5 0])
title([fig_name ' 3D projection plot'])

if params.objcolor <= 0.5
    colr = 'k';
elseif params.objcolor > 0.5
    colr = 'w';
end

%draw target outlines
outlines = cell(trajectory_frames,1);
outline_samples = 100;
for k = 1:trajectory_frames
    x_c = linspace(-a(k),a(k),outline_samples);
    y_c = sqrt(b(k).^2 - b(k).^2*x_c.^2 ./ a(k).^2);
    outline_2D = repmat(center(k,:),2*outline_samples,1)+[x_c fliplr(x_c) ; y_c -fliplr(y_c)]';
    ell = NaN(outline_samples,3);
    for i = 1:2*outline_samples
       ell(i,:) = feval(plane.getXYZ,outline_2D(i,1),outline_2D(i,2)); 
    end
    outlines{k} = ell;
end
plot3(ell(:,1),ell(:,2),ell(:,3),colr);
metrics.prey_outlines = outlines;


%% Find translation and rotation of the head

%calculate metric velocities and accelerations
[head_vel,head_accel,recon_head_vel,recon_head_accel,resid_head_vel,resid_head_accel,~,~]=numerical_differentiation(headPts,20,50,fps);
metrics.head_vel = head_vel; %convert to m*s^-1
metrics.head_accel = head_accel * 1000; %convert to m*s^-2

%Plot results
fig10 = figure;
subplot(3,1,1)
plot(time_vec,head_vel(:,1))
hold on
plot(time_vec,head_vel(:,2),'r')
plot(time_vec,head_vel(:,3),'g')
ylabel('Head velocity (m*s^{-1})')
legend('X','Y','Z')
title([fig_name ' Head velocity and acceleration'])

subplot(3,1,2)
plot(time_vec,head_accel(:,1))
hold on
plot(time_vec,head_accel(:,2),'r')
plot(time_vec,head_accel(:,3),'g')
ylabel('Head acceleration (m*s^{-2})')
legend('X','Y','Z')

subplot(3,1,3)
[ax,h1,h2] = plotyy(time_vec,sqrt(sum(head_vel.^2,2)),time_vec,sqrt(sum(head_accel.^2,2)));
set(get(ax(1),'Ylabel'),'String','Total head velocity (m*s^{-1})')
set(get(ax(2),'Ylabel'),'String','Total head acceleration (m*s^{-2})')
xlabel('Time (s)')
% set(h2,'Color','g')
saveas(fig10,[toad_folder filesep 'head velocity and acceleration.fig'])

%constuct head reference frame from head plane points
baseline = headPts(:,7:9) - headPts(:,4:6);
origin = repmat((dot(headPts(:,1:3) - headPts(:,7:9),baseline,2) ./ dot(baseline,baseline,2)),1,3).*baseline;
% x_h = headPts(:,4:6) - origin;
y_h = origin;
y_h = y_h ./ repmat(sqrt(sum(y_h.^2,2)),1,3);
x_h = headPts(:,1:3) - (origin + headPts(:,7:9));
x_h = x_h ./ repmat(sqrt(sum(x_h.^2,2)),1,3);
z_h = cross(x_h,y_h);
z_h = z_h ./ repmat(sqrt(sum(z_h.^2,2)),1,3);
metrics.head_frame = [x_h, y_h, z_h];

figure(fig1)
plot3([headPts(thit,1) headPts(thit,1)+x_h(thit,1)],[headPts(thit,2) headPts(thit,2)+x_h(thit,2)],[headPts(thit,3) headPts(thit,3)+x_h(thit,3)],'b')
plot3([headPts(thit,1) headPts(thit,1)+y_h(thit,1)],[headPts(thit,2) headPts(thit,2)+y_h(thit,2)],[headPts(thit,3) headPts(thit,3)+y_h(thit,3)],'r')
plot3([headPts(thit,1) headPts(thit,1)+z_h(thit,1)],[headPts(thit,2) headPts(thit,2)+z_h(thit,2)],[headPts(thit,3) headPts(thit,3)+z_h(thit,3)],'g')
saveas(fig1,[toad_folder filesep '3D tongue error.fig']);

%Calculate Z-Y-X Euler angles
phi = atan2d(x_h(:,2),x_h(:,1)) - atan2d(plane.normal(2),plane.normal(1)); %reference phi to screen plane normal
% phi = atan2d(x_h(:,2),x_h(:,1)) - atan2d(normal(2),normal(1)); %reference phi to screen plane normal
theta = -atan2d(x_h(:,3),sqrt(x_h(:,2).^2 + x_h(:,1).^2));
psi = atan2d(-y_h(:,3),-z_h(:,3));
metrics.euler_angles = [phi, theta, psi];

%Calculate angular velocity and acceleration of head
[euler_dt,euler_ddt,~,~,~,~,~,~] = numerical_differentiation([phi, theta, psi],20,50,fps);
omega_head(:,1) = euler_dt(:,2).*sind(phi) + euler_dt(:,3).*cosd(phi).*cosd(theta);
omega_head(:,2) = -euler_dt(:,2).*cosd(phi) + euler_dt(:,3).*sind(phi).*cosd(theta);
omega_head(:,3) = euler_dt(:,1) + euler_dt(:,3).*sind(theta);
metrics.head_angular_velocity = omega_head * 1000;

omega_head_dt(:,1) = euler_ddt(:,3) + sind(theta).*euler_ddt(:,1) + cosd(theta).*euler_dt(:,1).*euler_dt(:,2);
omega_head_dt(:,2) = cosd(theta).*sind(psi).*euler_ddt(:,1) - sind(theta).*sind(psi).*euler_dt(:,1).*euler_dt(:,2) + cosd(theta).*cosd(psi).*euler_dt(:,1).*euler_dt(:,3) - cosd(psi).*euler_ddt(:,2) + sind(psi).*euler_dt(:,2).*euler_dt(:,3);
omega_head_dt(:,3) = cosd(theta).*cosd(psi).*euler_ddt(:,1) - sind(theta).*cosd(psi).*euler_dt(:,1).*euler_dt(:,2) - cosd(theta).*sind(psi).*euler_dt(:,1).*euler_dt(:,3) + sind(psi).*euler_ddt(:,2) + cosd(psi).*euler_dt(:,2).*euler_dt(:,3);
metrics.head_angular_acceleration = omega_head_dt * 10^6;

% %calculate Euler angles
% x_p = cross(z_h,repmat([0 0 1],thit,1));
% alpha = atan2d(dot(x_p,y_h,2),dot(x_p,x_h,2));
% z_xy = hypot(z_h(:,1),z_h(:,2));
% beta = atan2d(z_xy,z_h(:,3));
% gamma = -atan2d(x_p(:,2),x_p(:,1));
%
% epsilon = 0.000001;
% idx = find(z_xy < epsilon);
% alpha(idx) = 0;
% beta(beta(idx) > 0) = 0;
% beta(beta(idx) <= 0) = pi;
% gamma(idx) = -atan2d(x_h(idx,2),x_h(idx,1));
%
% metrics.alpha = alpha;
% metrics.beta = beta;
% metrics.gamma = gamma;
%
% %calculate Tait-Bryan angles
% x_p = cross(z_h,repmat([0 1 0],thit,1));
% psi = acosd( dot(x_p,x_h,2)./ (sqrt(sum(x_p.^2,2))) );
% theta = 90 - acosd(z_h(:,2));
% phi = 90 - acosd(z_h(:,1));
%
% metrics.psi = psi;
% metrics.theta = theta;
% metrics.phi = phi;

%plot results
fig3 = figure;
subplot(3,1,1)
plot(time_vec,phi )
hold on
plot(time_vec,theta,'r')
plot(time_vec,psi,'g')
ylabel('Angular position (deg)')
legend('phi','theta','psi')
title([fig_name ' head angular velocity and acceleration'])

subplot(3,1,2)
plot(time_vec,omega_head(:,1))
hold on
plot(time_vec,omega_head(:,2),'r')
plot(time_vec,omega_head(:,3),'g')
ylabel('Angular velocity of head (deg*s^{-1})')
legend('X','Y','Z')

subplot(3,1,3)
plot(time_vec,omega_head_dt(:,1))
hold on
plot(time_vec,omega_head_dt(:,2),'r')
plot(time_vec,omega_head_dt(:,3),'g')
ylabel('Angular acceleration of head (deg*s^{-2})')
xlabel('Time (s)')
legend('X','Y','Z')
saveas(fig3,[toad_folder filesep 'Head angular velocity and accleration.fig']);

% fig3 = figure;
% subplot(3,1,1)
% plot(time_vec,alpha)
% hold on
% ylabel('\alpha (deg)')
% subplot(3,1,2)
% plot(time_vec,beta)
% hold on
% ylabel('\beta (deg)')
% subplot(3,1,3)
% plot(time_vec,gamma)
% xlabel('Time (sec)')
% ylabel('\gamma (deg)')
% saveas(fig3,[toad_folder filesep 'Euler angles plot.fig']);
%
% fig4 = figure;
% subplot(3,1,1)
% plot(time_vec,phi)
% hold on
% ylabel('\phi (deg)')
% subplot(3,1,2)
% plot(time_vec,theta)
% hold on
% ylabel('\theta (deg)')
% subplot(3,1,3)
% plot(time_vec,psi)
% xlabel('Time (sec)')
% ylabel('\psi (deg)')
% saveas(fig4,[toad_folder filesep 'Tait-Bryan angles plot.fig']);

%% Euclidian and angular error

threshold_window = 50;
norm_vel = sqrt(sum(head_vel(:,1:3).^2,2));
first_sample = find(~isnan(norm_vel));
if (first_sample(1)+threshold_window) > length(norm_vel)
    warning('Not enough samples of head movement. Aborting processing.');
    return
end
thresh = nanmean(norm_vel(first_sample(1):first_sample(1)+threshold_window))...
    +2*nanstd(norm_vel(first_sample(1):first_sample(1)+threshold_window));
moves = logical(norm_vel > thresh);
up_cross = find(moves(2:end) > moves(1:end-1));
if isempty(up_cross)
    up_cross = 1;
    %     'hold'
    %     return
end
t_0 = up_cross(end);
metrics.t_0 = t_0;
metrics.t_hit = thit;

for i = 1:3
    if i == 1
        c = metrics.target_coors.center(thit,:);
        f = 'centroid_error';
    elseif i == 2
        c = metrics.target_coors.leading_edge(thit,:);
        f = 'leading_edge_error';
    elseif i == 3
        c = metrics.target_coors.leading_edge(t_0,:);
        f = 't0_leading_edge_error';
    end
    error_vec = intersection-c;
    euclidian_error = norm(error_vec);
    %     t_u = (p_0 - p_1)./norm(p_0 - p_1);
    %     %     t_v = (p_0 - p_2)./norm(p_0 - p_2);
    %     if direction == 1
    %         t_v = cross(normal,t_u);
    %     elseif direction == -1
    %         t_v = cross(t_u,normal); %use normal to plane to ensure x and y axes are orthogonal
    %     end
    %     t_v = t_v./norm(t_v);
    %     x_err = t_u*((intersection-c)*t_u');
    %     y_err = t_v*((intersection-c)*t_v');
    %
    %     ts_err_x = t_u*((ts_intersection-c)*t_u');
    %     ts_err_y = t_v*((ts_intersection-c)*t_v');
    
    basis = plane.skew*plane.basis; %nb. skewed basis may not be orthogonal
    x_err = dot(error_vec,basis(:,1))./norm(basis(:,1));
    y_err = dot(error_vec,basis(:,2))./norm(basis(:,2));
    
    org = headPts(t_0,1:3);
    xz_tongue = (intersection-org) - y_h(t_0,:)*dot(intersection-org,y_h(t_0,:));
    xz_target = (c-org) - y_h(t_0,:)*dot(c-org,y_h(t_0,:));
    xy_tongue = (intersection-org) - z_h(t_0,:)*dot(intersection-org,y_h(t_0,:));
    xy_target = (c-org) - z_h(t_0,:)*dot(c-org,z_h(t_0,:));
    
    azimuthal_error = acosd(dot(xy_tongue,xy_target)/(norm(xy_tongue)*norm(xy_target)));
    elevational_error = acosd(dot(xz_tongue,xz_target)/(norm(xz_tongue)*norm(xz_target)));
    total_error = acosd(dot(intersection-org,c-org)/(norm(intersection-org)*norm(c-org)));
    
    metrics.euclidean_error.(f).x = x_err;
    metrics.euclidean_error.(f).y = y_err;
    metrics.euclidean_error.(f).total = euclidian_error;
    metrics.angular_error.(f).azimuth = azimuthal_error;
    metrics.angular_error.(f).elevation = elevational_error;
    metrics.angular_error.(f).elevation = total_error;
    
    %     xx = x_err + c;
    %     yy = y_err + c;
    %     tt = hypot(x_err,y_err) + c;
    %     tot_ang_err = acosd(dot(tt-headPts(t_0,1:3),c-headPts(t_0,1:3)) / (norm(tt-headPts(t_0,1:3))*norm(c-headPts(t_0,1:3))));
    %     x_ang_err = sign(x_err(2))*direction*acosd(dot(xx-headPts(t_0,1:3),c-headPts(t_0,1:3)) / (norm(xx-headPts(t_0,1:3))*norm(c-headPts(t_0,1:3))));
    %     y_ang_err = sign(-y_err(3))*acosd(dot(yy-headPts(t_0,1:3),c-headPts(t_0,1:3)) / (norm(yy-headPts(t_0,1:3))*norm(c-headPts(t_0,1:3))));
    
    %     metrics.(f) =  [euclidian_error, tot_ang_err;...
    %         norm(x_err)*sign(x_err(1))*-direction, x_ang_err;...
    %         norm(y_err)*sign(-y_err(3)), y_ang_err];
    %     metrics.(f)
end

% % Estimated parallax errors
% p_x = ts_err_x+c;
% p_y = ts_err_y+c;
% p_t = hypot(p_x,p_y)+c;
% parallax_total = acosd(dot(p_t-headPts(t_0,1:3),c-headPts(t_0,1:3)) / (norm(p_t-headPts(t_0,1:3))*norm(c-headPts(t_0,1:3))));
% parallax_x = metrics.centroid_error(2,2) - acosd(dot(p_x-headPts(t_0,1:3),c-headPts(t_0,1:3)) / (norm(p_x-headPts(t_0,1:3))*norm(c-headPts(t_0,1:3))));
% parallax_y = metrics.centroid_error(3,2) -  acosd(dot(p_y-headPts(t_0,1:3),c-headPts(t_0,1:3)) / (norm(p_y-headPts(t_0,1:3))*norm(c-headPts(t_0,1:3))));
% metrics.parallax_errors = [parallax_total; parallax_x; parallax_y];


%% Angular size and speed of target

pds = thit - length(metrics.target_coors.leading_edge);
if sign(pds) == -1
    pds = 0;
end

%angular size calculation
a = [NaN(pds,3); metrics.target_coors.leading_edge] - headPts(:,1:3);
b = [NaN(pds,3); metrics.target_coors.trailing_edge] - headPts(:,1:3);
c = [NaN(pds,3); metrics.target_coors.top_edge] - headPts(:,1:3);
d = [NaN(pds,3); metrics.target_coors.bottom_edge] - headPts(:,1:3);
metrics.prey_angular_size = [abs(acosd(dot(a,b,2)./ (sqrt(sum(a.^2,2)).*sqrt(sum(b.^2,2))) )),...
    abs(acosd(dot(c,d,2)./ (sqrt(sum(c.^2,2)).*sqrt(sum(d.^2,2))) ))]; %first position is size in azimuth, second elevation

%Target vector
c_vec = [NaN(pds,3); metrics.target_coors.center] - headPts(:,1:3);
norm_c_vec = sqrt(sum(c_vec.^2,2));
% % c_vec = [NaN(pds,3); target_coors] - origin;
metrics.prey_distance = norm_c_vec;
le_c_vec = [NaN(pds,3); metrics.target_coors.leading_edge] - headPts(:,1:3);

% %Rotate into toad frame using proper Euler angles
% c1 = cosd(metrics.alpha);
% c2 = cosd(metrics.beta);
% c3 = cosd(metrics.gamma);
% s1 = sind(metrics.alpha);
% s2 = sind(metrics.beta);
% s3 = sind(metrics.gamma);
% R = [c1(i)*c3(i) - c2(i)*s1(i)*s3(i), -c1(i)*s3(i) - c2(i)*c3(i)*s1(i), s1(i)*s2(i);...
%         c3(i)*s1(i) + c1(i)*c2(i)*s3(i), c1(i)*c2(i)*c3(i) - s1(i)*s3(i), -c1(i)*s2(i);...
%         s2(i)*s3(i), c3(i)*s2(i), c2(i)];
% trans_coors_eu = R*c_vec;
%Using Tait-Bryan angles

%project target vector onto head xy plane to find azimuth
proj_target = c_vec - repmat(dot(c_vec,z_h,2),1,3).*(z_h);
metrics.prey_azimuth = atan2d( dot(proj_target,y_h,2) , dot(proj_target,x_h,2) );
p_t = le_c_vec - repmat(dot(c_vec,z_h,2),1,3).*(z_h);
metrics.le_prey_azimuth = atan2d( dot(p_t,y_h,2) , dot(p_t,x_h,2) );
% metrics.prey_azimuth = acosd( dot(x_h,proj_target,2) ./ ( sqrt(sum(x_h.^2,2)) .* sqrt(sum(proj_target.^2,2)) ) );
% metris.prey_azimuth = atan2d( repmat(dot(y_h,proj_target,2),1,3).*y_h , repmat(dot(x_h,proj_target,2),1,3).*x_h );

%elevation by finding target vector angle w/ z axis
metrics.prey_elevation = -acosd( dot(z_h,c_vec,2) ./ ( sqrt(sum(z_h.^2,2)) .* norm_c_vec ) ) + 90;
metrics.le_prey_elevation = -acosd( dot(z_h,le_c_vec,2) ./ ( sqrt(sum(z_h.^2,2)) .* sqrt(sum(le_c_vec.^2,2)) ) ) + 90;

%Plot results
fig5 = figure;
[AX,H1,H2]=plotyy(time_vec,metrics.prey_azimuth,time_vec,metrics.prey_elevation);
hold on
set(get(AX(1),'Ylabel'),'String','Prey Azimuth (deg)')
set(get(AX(2),'Ylabel'),'String','Prey Elevation (deg)')
xlabel('Time (s)')
title([fig_name ' Azimuth and Elevation'])
line([0 time_vec(end)],[0 0],'Color','k','LineStyle','--')
% set(H2,'Color','r')
saveas(fig5,[toad_folder filesep 'Azimuth and elevation.fig'])

fig11 = figure;
[AX,H1,H2]=plotyy(time_vec,metrics.le_prey_azimuth,time_vec,metrics.le_prey_elevation);
hold on
set(get(AX(1),'Ylabel'),'String','Prey Azimuth (deg)')
set(get(AX(2),'Ylabel'),'String','Prey Elevation (deg)')
xlabel('Time (s)')
title([fig_name 'Leading Edge Azimuth and Elevation'])
line([0 time_vec(end)],[0 0],'Color','k','LineStyle','--')
% set(H2,'Color','r')
saveas(fig11,[toad_folder filesep 'Leading Edge Azimuth and elevation.fig'])

%prey angular velocity calculation
[prey_vel,prey_accel,prey_vel_recon,prey_accel_recon,prey_vel_resid,prey_accel_resid,~,~] = numerical_differentiation(metrics.target_coors.center,20,50,360);
prey_vel = prey_vel * 1000;
metrics.prey_theta = acosd(dot(prey_vel,c_vec,2) ./ ( sqrt(sum(prey_vel.^2,2)) .* norm_c_vec ) );
metrics.prey_angular_velocity = (180/pi) * cross(c_vec,[NaN(pds,3);prey_vel]) ./ repmat(norm_c_vec.^2,1,3);
metrics.prey_angular_speed = sqrt(sum(metrics.prey_angular_velocity.^2,2)).*sign(metrics.prey_angular_velocity(:,3)); %whether vector points up or down indicates left or rightward motion
norm_ang_vel = metrics.prey_angular_velocity ./ repmat(metrics.prey_angular_speed,1,3);
az_proj = norm_ang_vel - repmat(dot(norm_ang_vel,z_h,2),1,3).*(z_h);
% metrics.axis_of_rotation(:,1) = acosd( dot(x_h,az_proj,2) ./ sqrt(sum(az_proj.^2,2)) );
% metrics.axis_of_rotation(:,1) = atan2d( repmat(dot(y_h,az_proj,2),1,3).*y_h , repmat(dot(x_h,az_proj,2),1,3).*x_h );
metrics.axis_of_rotation(:,1) = atan2d( dot(y_h,az_proj,2) , dot(x_h,az_proj,2) );
metrics.axis_of_rotation(:,2) = -acosd(dot(norm_ang_vel,z_h,2)) + 90;  %defined relative to normal of head plane

%Angle of tongue and angular error
tongue_vec = intersection - headPts(t_0,1:3);
tongue_proj = tongue_vec - dot(tongue_vec,z_h(t_0,:)).*(z_h(t_0,:));
tongue_angle_az = atan2d( dot(tongue_proj,y_h(t_0,:)) , dot(tongue_proj,x_h(t_0,:)) );
tongue_angle_el =  -acosd( dot(z_h(t_0,:),tongue_vec) ./ norm(tongue_vec) ) + 90;
% metrics.tongue_error.angular.centroid.prey_coordinate_frame = ...
%     [(metrics.prey_azimuth(thit) - tongue_angle_az)*direction ,...
%     tongue_angle_el - metrics.prey_elevation(thit) ]; %errors projected onto planes defined by headPt(1) and prey trajectory
% metrics.tongue_error.angular.leading_edge.prey_coordinate_frame = ...
%     [(metrics.prey_azimuth(thit) - tongue_angle_az)*direction ,...
%     tongue_angle_el - metrics.prey_elevation(thit) ]; %errors projected onto planes defined by headPt(1) and prey trajectory
% metrics.tongue_error.angular.centroid.head_coordinate_frame = ...
%     [(metrics.prey_azimuth(thit) - tongue_angle_az)*direction ,...
%     tongue_angle_el - metrics.prey_elevation(thit) ]; %errors projected onto planes defined by head coordinate frame
% metrics.tongue_error.angular.leading_edge.head_coordinate_frame = ...
%     [(metrics.prey_azimuth(thit) - tongue_angle_az)*direction ,...
%     tongue_angle_el - metrics.prey_elevation(thit) ]; %errors projected onto planes defined by head coordinate frame

% metrics.centroid_angular_error_toad = [(metrics.prey_azimuth(thit) - tongue_angle_az)*direction , tongue_angle_el - metrics.prey_elevation(thit) ];
% metrics.leading_edge_angular_error_toad = [(metrics.le_prey_azimuth(thit) - tongue_angle_az)*direction , tongue_angle_el - metrics.le_prey_elevation(thit)  ];

%Plot results
fig6 = figure;
subplot(3,1,1)
plot(time_vec,metrics.prey_angular_size(:,1))
hold on
plot(time_vec,metrics.prey_angular_size(:,2),'r')
ylabel('Angular size (deg)')
xlabel('Time (s)')
legend('Azimuth','Elevation')
title([fig_name ' Prey angular velocity and acceleration'])
subplot(3,1,2)
[ax,h1,h2] = plotyy(time_vec,metrics.prey_angular_speed,time_vec,metrics.prey_distance);
set(get(ax(1),'Ylabel'),'String','Prey angular speed (deg/s)')
set(get(ax(2),'Ylabel'),'String','Prey distance (mm)')
% set(h2,'Color','r')
subplot(3,1,3)
plot(time_vec,metrics.axis_of_rotation(:,1))
hold on
plot(time_vec,metrics.axis_of_rotation(:,2),'r')
xlabel('Time (sec)')
ylabel('Axis of rotation angle (deg)')
legend('Azimuth','Elevation')
saveas(fig6,[toad_folder filesep 'Prey angular size and speed.fig']);

close all;

%% Metric size and speed of target

% metrics.prey_metric_size = params.objLenX * mmpp;
% % metrics.prey_metric_speed = params.objVelX * fps * mmpp;
% metrics.prey_metric_speed = sqrt(sum(prey_vel.^2,2));

% %% Subfunction 1
%     function[X,Y,Z] = warpPoints(points,plane)
%         
%         [x,y,z] = transformPointsForward(plane.rotate,points(:,1),points(:,2),points(:,3));
%         [v,w] = transformPointsForward(plane.skew,y,z);
%         [X,Y,Z] = transformPointsInverse(plane.rotate,x,v,w);
%         
%     end

end
