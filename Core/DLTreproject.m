function[uv] = DLTreproject(xyz,c)

%DLTREPROJECT Projects xyz points back onto images from multicam rig.
% Takes n x 3 array of xyz points, and 11 x m array of DLT coefficients
% from a m-camera rig. Returns n x 2 x m array of reprojected points.

nPts = size(xyz,1);
nCams = size(c,2);

uv = NaN(nPts,2,nCams);
for i = 1:nCams
    uv(:,1,i)=(xyz(:,1).*c(1,i)+xyz(:,2).*c(2,i)+xyz(:,3).*c(3,i)+c(4,i))./ ...
        (xyz(:,1).*c(9,i)+xyz(:,2).*c(10,i)+xyz(:,3).*c(11,i)+1);
    uv(:,2,i)=(xyz(:,1).*c(5,i)+xyz(:,2).*c(6,i)+xyz(:,3).*c(7,i)+c(8,i))./ ...
        (xyz(:,1).*c(9,i)+xyz(:,2).*c(10,i)+xyz(:,3).*c(11,i)+1);
end


