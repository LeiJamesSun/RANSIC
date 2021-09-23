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

function output_canvas = image_blending(warped_img1,warped_img2)

w1 = imfill(im2bw(uint8(warped_img1), 0),'holes');
w2 = imfill(im2bw(uint8(warped_img2), 0),'holes');

w1 = mat2gray(w1);
w2 = mat2gray(w2);

warped_img1 = double(warped_img1);
warped_img2 = double(warped_img2);
        
output_canvas(:,:,1) = ((warped_img1(:,:,1).*w1)+(warped_img2(:,:,1).*w2))./(w1+w2);
output_canvas(:,:,2) = ((warped_img1(:,:,2).*w1)+(warped_img2(:,:,2).*w2))./(w1+w2);
output_canvas(:,:,3) = ((warped_img1(:,:,3).*w1)+(warped_img2(:,:,3).*w2))./(w1+w2);

for i=1:5472
    for j=1:3648

        
        if max(output_canvas(i,j,:))~=0 && max(warped_img2(i,j,:))~=0
output_canvas(i,j,1) = ((warped_img2(i,j,1)));
output_canvas(i,j,2) = ((warped_img2(i,j,2)));
output_canvas(i,j,3) = ((warped_img2(i,j,3)));
        end
    end
end

% output_canvas(:,:,1) = ((warped_img1(:,:,1).*w1)+(warped_img2(:,:,1).*w2))./(w1+w2);
% output_canvas(:,:,2) = ((warped_img1(:,:,2).*w1)+(warped_img2(:,:,2).*w2))./(w1+w2);
% output_canvas(:,:,3) = ((warped_img1(:,:,3).*w1)+(warped_img2(:,:,3).*w2))./(w1+w2);
output_canvas = uint8(output_canvas);