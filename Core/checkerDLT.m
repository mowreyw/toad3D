function [plane] = checkerDLT(filename,camPts,checkerDim,origin,checkerSize)

% CHECKERDLT Get plane parameters for rear-projection screen from
% checkerboard image.
% 
% Takes filename of .csv file containing DLT coefficients for n cameras
% generated by DLTdv5. camPts is an p x 2n array of p checkerboard corners 
% (or other feature) observed in the n camera views. checkerDim describes the
% number of checkers in u and v. origin is the origin of the checkerboard
% in screen coordinates, and checkersize is the size of the checkers (or
% spacing of features) in uv coordinates.
% 
% Returns a structure 'plane' with the following fields:
% normal = normal vector for the plane
% basis = plane basis vectors aligned w/ u and v axes of the screen
% center_point = mean point of the featues in the projected image (lies on
% the plane)
% eig = eigenvectors of the feature cloud
% origin_uv = origin of the checkerboard in screen coordinates (u,v)
% origin_xyz = origin of the checkerboard in xyz coordinates (lies on the
% plane)
% DLT_coeffs = DLT coefficients
% xyz = xyz coordinates of the features
% rmse = root mean square reconstruction error
% uv = uv coordinates of the features
% 
% NB: when using 'calib_gui' to generate corners, be sure to correct the y
% coordinates by subtracting them from the image height.

% WRM 150102


% %constants for plotting the basis vectors
% org_pt = 10;
% norm_scale = 100;

%import DLT coeffs
c = importdata(filename);

%number of points
nPts = size(camPts,1);

% number of cameras
nCams=size(camPts,2)/2;

% setup output variables
xyz(1:nPts,1:3)=NaN;
rmse(1:nPts,1)=NaN;

% process each frame
for i=1:nPts

  % get a list of cameras with non-NaN [u,v]
  cdx=find(isnan(camPts(i,1:2:nCams*2))==false);

  % if we have 2+ cameras, begin reconstructing
  if numel(cdx)>=2

    % initialize least-square solution matrices
    m1=[];
    m2=[];

    m1(1:2:numel(cdx)*2,1)=camPts(i,cdx*2-1).*c(9,cdx)-c(1,cdx);
    m1(1:2:numel(cdx)*2,2)=camPts(i,cdx*2-1).*c(10,cdx)-c(2,cdx);
    m1(1:2:numel(cdx)*2,3)=camPts(i,cdx*2-1).*c(11,cdx)-c(3,cdx);
    m1(2:2:numel(cdx)*2,1)=camPts(i,cdx*2).*c(9,cdx)-c(5,cdx);
    m1(2:2:numel(cdx)*2,2)=camPts(i,cdx*2).*c(10,cdx)-c(6,cdx);
    m1(2:2:numel(cdx)*2,3)=camPts(i,cdx*2).*c(11,cdx)-c(7,cdx);

    m2(1:2:numel(cdx)*2,1)=c(4,cdx)-camPts(i,cdx*2-1);
    m2(2:2:numel(cdx)*2,1)=c(8,cdx)-camPts(i,cdx*2);

    % get the least squares solution to the reconstruction
    xyz(i,1:3)=linsolve(m1,m2);

    % compute ideal [u,v] for each camera
    uv=m1*xyz(i,1:3)';

    % compute the number of degrees of freedom in the reconstruction
    dof=numel(m2)-3;

    % estimate the root mean square reconstruction error
    rmse(i,1)=(sum((m2-uv).^2)/dof)^0.5;
  end
end

%fit plane to xyz points
[n,V,p] = affine_fit(xyz);
origin_mat = repmat(p,length(xyz),1);
xyz_in_plane = origin_mat + ...
    (xyz - origin_mat - repmat(n',length(xyz),1).*repmat((xyz-origin_mat)*n,1,3));
checker_points = reshape(xyz_in_plane,checkerDim(1),checkerDim(2),3);
direction_x = squeeze(mean(mean(diff(checker_points,1,1))));
norm_x = direction_x/norm(direction_x);
% x_in_plane = direction_x - n*(direction_x'*n);
% norm_x = x_in_plane/norm(x_in_plane);
% norm_y = cross(-n,norm_x);
norm_y = cross(norm_x,n);
% if sign(n(1))
%     norm_y = cross(n,norm_x);
% else
%     norm_y = cross(norm_x,n);  
% end
% norm_y = cross(n,norm_x);
% direction_y = squeeze(mean(mean(diff(checker_points,1,2))));
% basis = [direction_x/norm(direction_x) direction_y/norm(direction_y)]; %should enforce orthogonality here
basis = [norm_x norm_y];
% dist = (squeeze(checker_points(1,1,:))' - p)*n;

%save data
plane.normal = n;
plane.basis = basis;
plane.eig = V;
plane.center_point = p;
plane.origin_uv = origin;
plane.origin_xyz = xyz_in_plane(prod(checkerDim) - (checkerDim(1) -1),:); 
plane.DLT_coeffs = c;
plane.xyz = xyz;
plane.rmse = rmse;
% plane.mmpp = norm(xyz(checkerDim(1),:) -  xyz(1,:)) / (checkerSize*(checkerDim(1)-1));
plane.mmpp = [norm(xyz(checkerDim(1),:) -  xyz(1,:)) / (checkerSize*(checkerDim(1)-1))...
     norm(xyz(prod(checkerDim),:) -  xyz(checkerDim(1)-1,:)) / (checkerSize*(checkerDim(2)-1))...
     1];

%build matrices for moving between uv and xyz coordinates
uv2xyz = [basis plane.origin_xyz'];
% xyz2uv = inv(uv2xyz);

%generate xyz coordinates from transformation matrix
checkU = repmat(origin(1)+(0:checkerSize:(checkerDim(1)-1)*checkerSize)',checkerDim(2),1);
checkV = reshape(repmat(origin(2)+(0:checkerSize:(checkerDim(2)-1)*checkerSize),checkerDim(1),1),...
    prod(checkerDim),1);
checkV = flipud(checkV); %image points are ordered from top left corner down
recon_xyz = NaN(length(checkU),3);
for i = 1:length(checkU)
    recon_xyz(i,:) = uv2xyz*(([checkU(i) checkV(i) 1]'-[plane.origin_uv 0]').*plane.mmpp');  
end

%direct linear tranformation to correct skew
source_points = NaN(3*length(recon_xyz),9);
for i = 1:length(recon_xyz)
    source_points((i-1)*3+1:i*3,:) = blkdiag(recon_xyz(i,:),recon_xyz(i,:),recon_xyz(i,:));
end
a =  source_points \ reshape(xyz_in_plane',3*length(recon_xyz),1);
A = reshape(a,3,3)';
warped = NaN(length(recon_xyz),3);
for i = 1:length(recon_xyz)
    warped(i,:) = A*recon_xyz(i,:)';
end

plane.uv2xyz = uv2xyz;
% plane.xyz2uv = xyz2uv; %save basis and skew matrices separately to see errors in unskewed data
plane.skew = A;
plane.uv2xyz_skew = A*uv2xyz;
% plane.getUV = @(q) ((xyz2uv*q)./plane.mmpp') + [plane.origin_uv 0]';
plane.getUV = @(q) ((plane.uv2xyz_skew\q)./plane.mmpp') + [plane.origin_uv 0]';
plane.getXYZ = @(u,v) plane.uv2xyz_skew*(([u v 1]'-[plane.origin_uv 0]').*plane.mmpp');
save([pwd filesep 'plane.mat'],'plane');

figure;
scatter3(xyz_in_plane(:,1),xyz_in_plane(:,2),xyz_in_plane(:,3),'bo')
hold on
scatter3(recon_xyz(:,1),recon_xyz(:,2),recon_xyz(:,3),'ro')
scatter3(warped(:,1),warped(:,2),warped(:,3),'go')
legend('Measured corners','Corners from plane eq','Warp to fit')

end




% %correct for skew in uv coordinates
% K1 = [0 0 -1; 0 0 0; 1 0 0]; %rotation around y axis to align normal with x
% K2 = [0 -1 0; 1 0 0; 0 0 0]; %rotation around z to align normal with x
% theta = acos(dot([1 0],[n(1) n(3)])/norm([n(1) n(3)])); %calculate angle between normal and x axis in xz plane
% phi = acos(dot([1 0],n(1:2))/norm(n(1:2))); %calculate angle between normal and x in xy plane
% rot = [[expm(theta*K1)*expm(phi*K2) -p']; 0 0 0 1];
% inv_rot = inv(rot);
% checkU = repmat(origin(1)+(0:checkerSize:(checkerDim(1)-1)*checkerSize)',checkerDim(2),1);
% checkV = reshape(repmat(origin(2)+(0:checkerSize:(checkerDim(2)-1)*checkerSize),checkerDim(1),1),...
%     prod(checkerDim),1);
% checkV = flipud(checkV); %image points are ordered from top left corner down
% 
% recon_xyz = NaN(length(checkU),3);
% data_yz = NaN(length(checkU),4);
% recon_yz =data_yz;
% for i = 1:length(checkU)
%     recon_xyz(i,:) = uv2xyz*(([checkU(i) checkV(i) 1]'-[plane.origin_uv 0]').*plane.mmpp');
%     data_yz(i,:) = (rot*[xyz_in_plane(i,:) 1]')'; 
%     recon_yz(i,:) = (rot*[recon_xyz(i,:) 1]')';    
% end
% [Y,Z] = estimateGeometricTransform();
% [Y3,Z3] = transformPointsForward(tform,Y2,Z2);
% [X4,Y4,Z4] = transformPointsInverse(trans,X2,Y3,Z3);




% %save data
% plane.normal = n;
% plane.basis = basis;
% plane.eig = V;
% plane.center_point = p;
% plane.origin_uv = origin;
% plane.origin_xyz = xyz_in_plane(1,:); 
% plane.DLT_coeffs = c;
% plane.xyz = xyz;
% plane.rmse = rmse;
% % plane.mmpp = norm(xyz(checkerDim(1),:) -  xyz(1,:)) / (checkerSize*(checkerDim(1)-1));
% plane.mmpp = [norm(xyz(checkerDim(1),:) -  xyz(1,:)) / (checkerSize*(checkerDim(1)-1))...
%      norm(xyz(prod(checkerDim),:) -  xyz(checkerDim(1)-1,:)) / (checkerSize*(checkerDim(2)-1))...
%      1];
% % plane.eq = @(u,v) repmat(plane.origin_xyz,length(u),1) + ...
% %     plane.mmpp*repmat(u-plane.origin_uv(1),1,3).*repmat(plane.basis(:,1)',length(u),1)...
% %     + plane.mmpp*repmat(v-plane.origin_uv(2),1,3).*repmat(plane.basis(:,2)',length(u),1);
% % plane.getUV = @(q) xyz2uv*((q-plane.origin_xyz')./plane.mmpp');
% plane.getUV = @(q) ((xyz2uv*q)./plane.mmpp') + [plane.origin_uv 0]';
% plane.getXYZ = @(u,v) uv2xyz*(([u v 1]'-[plane.origin_uv 0]').*plane.mmpp');
% 
% %get uv points of checkerboard
% checkU = repmat(origin(1)+(0:checkerSize:(checkerDim(1)-1)*checkerSize)',checkerDim(2),1);
% checkV = reshape(repmat(origin(2)+(0:checkerSize:(checkerDim(2)-1)*checkerSize),checkerDim(1),1),...
%     prod(checkerDim),1);
% checkV = flipud(checkV); %image points are ordered from top left corner down
% % end_points = origin + checkerSize.*(checkerDim-1);
% % checkX = repmat(linspace(origin(1),end_points(1),checkerDim(1))',checkerDim(2),1);
% % [checkX,checkY] = meshgrid(linspace(origin(1),end_points(1),checkerDim(1)),...
% %     linspace(origin(2),end_points(2),checkerDim(2)));
% % checkX = flipud(checkX');
% % checkY = flipud(checkY');

% %test of transformation
% for i = 1:length(checkU)
%     recon_xyz(i,:) = feval(plane.getXYZ,checkU(i),checkV(i))';
%     uv_reproj(i,:) = feval(plane.getUV,recon_xyz(i,:)');
% end

% % compare uv to reconstructed uv and find transformations to correct skew
% recon_xyz = feval(plane.eq,checkU(:),checkV(:));
% K1 = [0 0 1; 0 0 0; -1 0 0]; %rotation around y axis to align normal with x
% K2 = [0 -1 0; 1 0 0; 0 0 0]; %rotation around z to align normal with x
% theta = acos(dot([1 0],[n(1) n(3)])); %calculate angle between normal and x axis in xz plane
% phi = acos(dot([1 0],n(1:2))); %calculate angle between normal and x in xy plane
% rot = expm(theta*K1)*expm(phi*K2);
% trans = affine3d([[rot zeros(3,1)]; [-p 1]]);
% [X1,Y1,Z1] = transformPointsForward(trans,xyz_in_plane(:,1),xyz_in_plane(:,2),xyz_in_plane(:,3));
% [X2,Y2,Z2] = transformPointsForward(trans,recon_xyz(:,1),recon_xyz(:,2),recon_xyz(:,3));
% tform = estimateGeometricTransform([Y2 Z2],[Y1 Z1],'projective');
% [Y3,Z3] = transformPointsForward(tform,Y2,Z2);
% [X4,Y4,Z4] = transformPointsInverse(trans,X2,Y3,Z3);
% % y_rot = affine3d([[expm(theta*K1) zeros(3,1)];[0 0 0 1]]);
% % [X1,Y1,Z1] = transformPointsForward(y_rot,xyz_in_plane(:,1)-p(1),xyz_in_plane(:,2)-p(2),xyz_in_plane(:,3)-p(3));
% % y_rot = affine3d([[expm(theta*K1) zeros(3,1)];[-p 1]]);
% % [X1,Y1,Z1] = transformPointsForward(y_rot,xyz_in_plane(:,1),xyz_in_plane(:,2),xyz_in_plane(:,3));
% % [X2,Y2,Z2] = transformPointsForward(y_rot,recon_xyz(:,1),recon_xyz(:,2),recon_xyz(:,3));
% % tform = estimateGeometricTransform([Y2 Z2],[Y1 Z1],'projective');
% % [Y3,Z3] = transformPointsForward(tform,Y2,Z2);
% % [X4,Y4,Z4] = transformPointsInverse(y_rot,X2,Y3,Z3);
% % [X4,Y4,Z4] = transformPointsInverse(y_rot,X2,Y3,Z3);

% figure;
% scatter3(xyz_in_plane(:,1),xyz_in_plane(:,2),xyz_in_plane(:,3),'bo')
% hold on
% scatter3(recon_xyz(:,1),recon_xyz(:,2),recon_xyz(:,3),'ro')
% scatter3(X4,Y4,Z4,'go')
% legend('Measured corners','Corners from plane eq','Warp to fit')

% % save data
% plane.rotate = y_rot;
% plane.skew = tform;
% 
% save([pwd filesep 'plane.mat'],'plane');


% %find uv coordinates of xyz corners
% shifted_plane = xyz - repmat(xyz(1,:),length(xyz),1);
% K1 = [0 0 1; 0 0 0; -1 0 0]; %rotation around y axis to align normal with reverse x 
% theta = acos(dot([1 0 0],n)); %calculate angle between normal and reverse x axis
% y_rot = affine3d([[expm(theta*K1) zeros(3,1)];[0 0 0 1]]);
% [X1,Y1,Z1] = transformPointsForward(y_rot,shifted_plane(:,1),...
%     shifted_plane(:,2),shifted_plane(:,3));
% % K2 = [0 0 0; 0 0 1; 0 -1 0]; %rotation about reverse x axis to align norm_x with y
% cp = reshape([Y1 Z1],checkerDim(1),checkerDim(2),2);
% % cp = reshape([X1 Y1 Z1],checkerDim(1),checkerDim(2),3);
% direction_u = squeeze(mean(mean(diff(cp,1,1))));
% norm_u = direction_u/norm(direction_u);
% phi = acos(dot([1 0],norm_u)); %calculate angle between u and y
% x_rot = affine2d([cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]);
% [Y2,Z2] = transformPointsForward(x_rot,Y1,Z1);
% exp(theta+phi
% % phi = acos(dot([0 1 0],[0 norm_u(2) norm_u(3)]));
% phi = acos(dot([0 1 0],norm_x)); 
% % x_rot = affine3d([[expm(phi*K2) zeros(3,1)];[0 0 0 1]]);
% % [X2,Y2,Z2] = transformPointsForward(x_rot,X1,Y1,Z1); %is there a way to compose these transforms into one?
% u_scale = (checkerDim(1)*checkerSize)/Y2(checkerDim(1));
% % v_scale = rotated_plane(1+checkerDim(1)*(checkerDim(2)-1),2);
% scale = affine2d([u_scale*eye(2) [0; 0]; 0 0 1]);
% scaled_plane = repmat(origin,prod(checkerDim),1) + transformPointsForward(scale,[Y2 Z2]);
% 
% % find projective transformation from regenerated corners to actual corners
% tform = estimateGeometricTransform([checkX(:) checkY(:)],scaled_plane,'projective');
% [Y3,Z3] = transformPointsInverse(tform,scaled_plane(:,1),scaled_plane(:,2));

%

% xyz = [xyz; xyz(org_pt,:) + norm_scale*norm_x'; xyz(org_pt,:) + norm_scale*norm_y'];
% xyz = [xyz; xyz(org_pt,:) + norm_scale*V(:,1)'; xyz(org_pt,:) + norm_scale*V(:,2)'];

% %check quality of solution by reprojecting points onto original images
% cam1_uv(:,1)=(xyz(:,1).*c(1,1)+xyz(:,2).*c(2,1)+xyz(:,3).*c(3,1)+c(4,1))./ ...
%   (xyz(:,1).*c(9,1)+xyz(:,2).*c(10,1)+xyz(:,3).*c(11,1)+1);
% cam1_uv(:,2)=(xyz(:,1).*c(5,1)+xyz(:,2).*c(6,1)+xyz(:,3).*c(7,1)+c(8,1))./ ...
%   (xyz(:,1).*c(9,1)+xyz(:,2).*c(10,1)+xyz(:,3).*c(11,1)+1);
% cam2_uv(:,1)=(xyz(:,1).*c(1,2)+xyz(:,2).*c(2,2)+xyz(:,3).*c(3,2)+c(4,2))./ ...
%   (xyz(:,1).*c(9,2)+xyz(:,2).*c(10,2)+xyz(:,3).*c(11,2)+1);
% cam2_uv(:,2)=(xyz(:,1).*c(5,2)+xyz(:,2).*c(6,2)+xyz(:,3).*c(7,2)+c(8,2))./ ...
%   (xyz(:,1).*c(9,2)+xyz(:,2).*c(10,2)+xyz(:,3).*c(11,2)+1);
% 
% %plot to check the basis vectors for alignment with u and v axes
% figure
% scatter3(xyz(1:60,1),xyz(1:60,2),xyz(1:60,3),'r.')
% hold on
% plot3([xyz(org_pt,1) xyz(end-3,1)],[xyz(org_pt,2) xyz(end-3,2)],[xyz(org_pt,3) xyz(end-3,3)],'b')
% plot3([xyz(org_pt,1) xyz(end-2,1)],[xyz(org_pt,2) xyz(end-2,2)],[xyz(org_pt,3) xyz(end-2,3)],'g')
% plot3([xyz(org_pt,1) xyz(end-1,1)],[xyz(org_pt,2) xyz(end-1,2)],[xyz(org_pt,3) xyz(end-1,3)],'b')
% plot3([xyz(org_pt,1) xyz(end,1)],[xyz(org_pt,2) xyz(end,2)],[xyz(org_pt,3) xyz(end,3)],'g')
% legend('corners','norm_x','norm_y','V_x','V_y')
% axis equal


