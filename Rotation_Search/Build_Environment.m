function [pts_3d,pts_3d_,R]=Build_Environment(n_ele,noise,outlier_ratio)

pts_3d=rand(n_ele,3)-0.5;

%% Transformation

% rotation
axis1= rand(1,3)-0.5;

axis1=axis1/norm(axis1);

ag=2*pi*rand(1);

A1 = ag * axis1;

if norm(A1) < 2e-16
         R=eye(3);
else a = A1 / norm(A1);
K1=[0, -a(3), a(2); ...
    a(3), 0, -a(1); ...
    -a(2), a(1), 0];

R = eye(3) + sin(norm(A1)) * K1 + (1 -cos(norm(A1))) * K1 * K1;
end


for i=1:n_ele
    pts_3d(i,:)=pts_3d(i,:)/norm(pts_3d(i,:));
end

%transform by R & t
pc_med= pts_3d * R';

% add noise
    for i=1:n_ele
    pc_med(i,:)=pc_med(i,:)+noise*randn(1,3);
    pts_(i,:) = pc_med(i,:)/norm(pc_med(i,:));
    end

pts_3d_=pts_;

% create outliers
for i=1:round(n_ele*outlier_ratio)
    
    for iii=1:1e+18
    rand_vec=2*rand(1,3)-1;
        if norm(rand_vec)<=0.5*sqrt(3)
            break
        end
    end
    
pts_3d_(i,:)=rand_vec/norm(rand_vec);

end



%% Show Figure

show_figure=0;

if show_figure==1
    
    for i=1:n_ele
    plot3([0,pts_3d(i,1)],[0,pts_3d(i,2)],[0,pts_3d(i,3)],'r','LineWidth',1);
    hold on;
    end
    
end


% if show_figure==1
%     
% figure;
% 
% pc1=pointCloud(pts_3d(1:n_ele,:));
% pc2_out=pointCloud(pts_3d_(1:round(outlier_ratio*n_ele),:));
% pc2_in=pointCloud(pts_3d_(1+round(outlier_ratio*n_ele):n_ele,:));
% 
% pcshow([pc1.Location(:,1),pc1.Location(:,2),pc1.Location(:,3)],[0 0 1],'MarkerSize',90);
% 
% hold on;
% 
% pcshow([pc2_out.Location(:,1),pc2_out.Location(:,2),pc2_out.Location(:,3)],[1 0 0],'MarkerSize',200);
% 
% hold on;
% 
% pcshow([pc2_in.Location(:,1),pc2_in.Location(:,2),pc2_in.Location(:,3)],[0 1 0],'MarkerSize',200);
% 
% hold on;
% 
% for i=round(n_ele*outlier_ratio)+1:n_ele
%     
%     plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'g','LineWidth',2.5);
%     
% end
% 
% for i=1:round(n_ele*outlier_ratio)
%     
%     pp=plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'r','LineWidth',1);
% 
%     pp.Color(4)=0.15;
%     
% end
% 
% grid off;
% axis off;
% 
% end

end