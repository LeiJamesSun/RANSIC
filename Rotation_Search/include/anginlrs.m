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


function [ idx ] = anginlrs(X,Y,th)
D = normr(X)-normr(Y);
len = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
idx = len<=th; %2*sin(th/2);
end

