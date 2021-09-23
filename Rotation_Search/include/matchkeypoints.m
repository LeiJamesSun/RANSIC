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

function [ match_M, match_B, D, KPD, idx_M, idx_B ] = ...
    matchkeypoints(M, B, desc_M, desc_B, th)
% The result is sorted 

msize=length(M);
bsize=length(B);

match_M=zeros(0,3);
match_B=zeros(0,3);
idx_M=[];
idx_B=[];
D=[];
KPD=[];

for j=1:bsize
    b = B(j, :);
    desc_b = desc_B(j, :);
    
    for i=1:msize
            m = M(i, :);
            desc_m = desc_M(i, :);
            
            ddist = norm(desc_m-desc_b);
            %pdist = norm(m-b);
        
        if ddist < th
            match_M(end+1,:)=m;
            match_B(end+1,:)=b;
            idx_M(end+1)=i;
            idx_B(end+1)=j;
            D(end+1)=ddist;
            KPD(end+1)=norm(m-b);
        end
    end
end
[D, idx]= sort(D);
match_M = match_M(idx,:);
match_B = match_B(idx,:);
idx_M = idx_M(idx);
idx_B = idx_B(idx);
KPD = KPD(idx);

end

