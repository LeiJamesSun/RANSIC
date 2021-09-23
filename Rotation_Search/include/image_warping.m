%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     G U A R A N T E E D    O U T L I E R    R E M O V A L
%            F O R   R O T A T I O N   S E A R C H
%
%
% This package contains the source code which implements the
% guaranteed outlier removal for rotation search proposed in
% Alvaro PARRA BUSTOS, Tat-Jun CHIN
% Guaranteed Outlier Removal for Rotation Search
% In International Conference on Computer Vision (ICCV) Dec 2015, Santiago
%
% Copyright (c) 2015 Alvaro PARRA BUSTOS (aparra@cs.adelaide.edu.au.)
% School of Computer Science, The University of Adelaide, Australia
% The Australian Center for Visual Technologies
% http://cs.adelaide.edu.au/~aparra
% Please acknowledge the authors by citing the above paper in any academic
% publications that have made use of this package or part of it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ imgwarped ] = image_warping(img1,img2,H12)

imgg2 = rgb2gray(uint8(img2));


[height1,width1,~] = size(img1);
[height2,width2,~] = size(img2);

imgwarped = zeros(size(img1));

[~,~,d] = size(img1);
[posx,posy] = meshgrid(1:width1,1:height1);
pos1 = [posx(:)';posy(:)'];
pos2 = H12 * makehomogeneous(pos1);
pos2 = round(makeinhomogeneous(pos2));
in = inpolygon(pos2(1,:),pos2(2,:),[1 width1 width1 1 1],[1 1 height1 height1 1]);
pos1 = pos1(:,in);
pos2 = pos2(:,in);

in = inpolygon(pos2(1,:),pos2(2,:),[1 width2 width2 1 1],[1 1 height2 height2 1]);
pos1 = pos1(:,in);
pos2 = pos2(:,in);


img2 = double(img2);
[y,x] = find(imgg2 == 0);

if d==3
    for i = 1:length(x)
        img2(y(i),x(i),1) = 0;
        img2(y(i),x(i),2) = 0;
        img2(y(i),x(i),3) = 0;
    end
elseif d==1
    for i = 1:length(x)
        img2(y(i),x(i)) = 0;
        img2(y(i),x(i)) = 0;
        img2(y(i),x(i)) = 0;
    end
end

if d==3
    for i = 1:size(pos1,2);
        imgwarped(pos1(2,i),pos1(1,i),1) = img2(pos2(2,i),pos2(1,i),1);
        imgwarped(pos1(2,i),pos1(1,i),2) = img2(pos2(2,i),pos2(1,i),2);
        imgwarped(pos1(2,i),pos1(1,i),3) = img2(pos2(2,i),pos2(1,i),3);
    end
elseif d==1
    for i = 1:size(pos1,2);
        imgwarped(pos1(2,i),pos1(1,i)) = img2(pos2(2,i),pos2(1,i));
        imgwarped(pos1(2,i),pos1(1,i)) = img2(pos2(2,i),pos2(1,i));
        imgwarped(pos1(2,i),pos1(1,i)) = img2(pos2(2,i),pos2(1,i));
    end
end
if d == 3
    imgwarped = uint8(imgwarped);
end

end
