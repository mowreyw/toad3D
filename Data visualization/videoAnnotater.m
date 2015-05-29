function[]=videoAnnotater(directory,c)

%VIDEOANNOTATER Annotate videos with 3D tracking results.

marker_size = 4;
fig_width = 2000;
fig_height = 1000;
trailing_frames = 20; %hold on tongue projection at end

%% Find index of toad movies
if ~strcmp(pwd,directory)
    cd(directory);
end

if strcmp(computer,'MACI64')
    toadDir = ls('-d',[pwd filesep '*b*cam2*']);
else
    toadDir = ls([pwd filesep '*b*cam2*']);
end

%% annotate

for k = 1:size(toadDir,1)
    
    filename = deblank(toadDir(k,:));
    writer = VideoWriter([pwd filesep filename filesep filename ' annotated video.avi']);
    open(writer);
    if ~exist([filename filesep 'metrics.mat'],'file')
        continue
    else
        load([filename filesep 'metrics.mat']);
    end
    target.center = DLTreproject(metrics.target_coors.center,c);
    target.center_raw = DLTreproject(metrics.target_coors_raw.center,c);
%     target.leading_edge = DLTreproject(metrics.target_coors.leading_edge,c);
%     target.top = DLTreproject(metrics.target_coors.top_edge,c);
    toad.pt1 = DLTreproject(metrics.headPts(:,1:3),c);
    toad.pt2 = DLTreproject(metrics.headPts(:,4:6),c);
    toad.pt3 = DLTreproject(metrics.headPts(:,7:9),c);
    
    nCams = size(target.center,3);
    
    %loop to write annotated video
    f1 = figure;
    set(f1,'Position',[100,100,fig_width,fig_height])
    vid_title = strrep(filename,'_',' ');
    for i = 1:metrics.t_hit
        po = DLTreproject(metrics.prey_outlines{i},c);
        for j = 1:nCams
            idx = num2str(j);
            if i == 1
                cam_path = strrep(filename,'cam2',['cam' idx]);
                readers.(['cam' idx]) = VideoReader([cam_path filesep cam_path '.avi']);
            end
            im = read(readers.(['cam' idx]),i);
            subplot(1,nCams,j)
            hold off
%             imshow(im)
            imshow(imadjust(im,[],[],0.7));
            hold on
            plot(target.center(i,1,j),readers.(['cam' idx]).Height - target.center(i,2,j),'yo','MarkerSize',marker_size)
            plot(target.center(1:i,1,j),readers.(['cam' idx]).Height - target.center(1:i,2,j),'y')
            plot(po(:,1,j),readers.(['cam' idx]).Height - po(:,2,j),'r')
%             plot(target.center_raw(i,1,j),readers.(['cam' idx]).Height - target.center_raw(i,2,j),'mo','MarkerSize',marker_size)
%             plot(target.center_raw(1:i,1,j),readers.(['cam' idx]).Height - target.center_raw(1:i,2,j),'m')
%             plot(target.leading_edge(i,1,j),readers.(['cam' idx]).Height - target.leading_edge(i,2,j),'co','MarkerSize',marker_size)
%             plot(target.top(i,1,j),readers.(['cam' idx]).Height - target.top(i,2,j),'mo','MarkerSize',marker_size)
            plot(toad.pt1(i,1,j),readers.(['cam' idx]).Height - toad.pt1(i,2,j),'bo','MarkerSize',marker_size)
            plot(toad.pt2(i,1,j),readers.(['cam' idx]).Height - toad.pt2(i,2,j),'ro','MarkerSize',marker_size)
            plot(toad.pt3(i,1,j),readers.(['cam' idx]).Height - toad.pt3(i,2,j),'go','MarkerSize',marker_size)
            plot(toad.pt1(1:i,1,j),readers.(['cam' idx]).Height - toad.pt1(1:i,2,j),'b')
            plot(toad.pt2(1:i,1,j),readers.(['cam' idx]).Height - toad.pt2(1:i,2,j),'r')
            plot(toad.pt3(1:i,1,j),readers.(['cam' idx]).Height - toad.pt3(1:i,2,j),'g')
            if i == metrics.t_hit
                %tongue projection code
            end
            title([vid_title ' tracking video, frame ' num2str(i)])
        end
        frame = getframe(f1);
        writeVideo(writer,frame);
    end
    
    for l = 1:trailing_frames
        writeVideo(writer,frame);
    end
    close(writer);
    close(f1);
    
end