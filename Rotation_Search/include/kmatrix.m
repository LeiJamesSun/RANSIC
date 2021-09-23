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


function [K]=kmatrix(filename, I)

%focal length in pixels = 
% (image width in pixels) * (focal length in mm) / (CCD width in mm)

models   = {'DSC-W1', 'Canon PowerShot A540','DSC-P100' }; 
values = [ 7.144, 5.744, 7.144];
Map = containers.Map(models,values);

info = imfinfo(filename);
model=info.Model;

% Calibration matrix 1
w=info.DigitalCamera.CPixelXDimension;
fl=info.DigitalCamera.FocalLength;

if ~isKey(Map, model)
    error('Add CCD width for camera model %s\n%s\n', model, ...
    sprintf('http://www.dpreview.com/search?query="%s"', model));
    
end

ccdw = Map(model);

f=w*fl/ccdw;

[r,c,~]=size(I);
K=zeros(3,3);
K(1,1)=f;
K(2,2)=f;
K(3,3)=1;
K(1,3)=c/2;
K(2,3)=r/2;

end