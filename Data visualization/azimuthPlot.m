function[] = azimuthPlot(varargin)

%TOAD3D Takes directory of toad data files and saves tongue error data and
%figures

f1=figure;
f2=figure;

prey_az = [];
prey_spd = [];
for j = 1:nargin
    
    toadDir = dir(varargin{j});
    
    for i = 1:size(toadDir,1)
        toad_folder = [varargin{j} filesep deblank(toadDir(i).name)];
        fprintf('%s\n',['Processing ' toad_folder]);
        if ~exist([toad_folder filesep 'metrics.mat'],'file')
            fprintf('%s\n','No metrics file found');
            continue
        end
        
        load([toad_folder filesep 'metrics.mat'])
        load(deblank(ls('-d',[toad_folder filesep 'stim param*'])))
        if ~isfield(metrics,'t_0')
            continue
        end
        
        c = colormap(jet);
%         stim_vec = 64.*[40:40:200]./200;
        time_vec = [-metrics.t_0:1:0]./metrics.fps;
        idx = find(~isnan(metrics.le_prey_azimuth));
        if sign(metrics.prey_azimuth(idx(1))) == -1
            az_data = -metrics.le_prey_azimuth;
        else
            az_data = metrics.le_prey_azimuth;
        end
        
%         plot(time_vec,az_data(1:metrics.t_0+1),'Color',c(round(64*params.objLenX/200),:))
%         plot(time_vec,metrics.prey_azimuth(1:metrics.t_0+1))
        figure(f1);
        plot(time_vec,az_data(1:metrics.t_0+1),'Color',c(round(64*metrics.prey_angular_speed(metrics.t_0)./100),:))
        hold on
        
        prey_spd = [prey_spd; metrics.prey_angular_speed(metrics.t_0)];
        prey_az = [prey_az; abs(metrics.le_prey_azimuth(metrics.t_0))];
       
    end
    
end
xlabel('Time (s)')
ylabel('Azimuth (deg)')

figure(f2);
scatter(prey_spd,prey_az,'filled')
xlabel('Prey angular speed (deg/s)')
ylabel('Prey azimuth (deg0')