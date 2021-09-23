%% Solver 'RANSIC' for Point Cloud Registration

% This is the source code of paper: 'RANSIC: Fast and Highly Robust 
% Estimation for Rotation Search and Point Cloud Registration using 
% Invariant Compatibility'

% Copyright by Lei Sun (leisunjames@126.com)

clc;
clear all;
close all;

%% Parameter Setup by User(s): 

% whether the scale is known
known_scale=0; % '0' for unknown scale and '1' for known scale (s=1)

% set outlier ratio in percentage
outlier_ratio=0.90; % should be set within 0-0.99 (0.99 could be slower)

% whether to display the point correspondences
show_figure=1; % '1' for displaying and '0' for not displaying



%% RANSIC Begins: (No Need for Adjustment)

disp(['RANSIC for Point Cloud Registration :']);

% set point correspondence number
n_ele=1000;

disp(['N=',num2str(n_ele)]);

% set noise in pixels
noise=0.01;

disp(['Noise=',num2str(noise)]);

% scale is within 1-5 (unknown scale) or 1 (known scale)
if known_scale==1 %known scale
scale_gt=1;
elseif known_scale==0 % unknown scale
scale_gt=1+4*rand(1);
end

% build simulation environment

% Input:
% pts_3d is a Nx3 matrix
% pts_3d_ is a Nx3 matrix
% R_gt is a 3x3 SO(3) matrix ('gt' means ground-truth)
% t_gt is a 1x3 vector

[pts_3d,pts_3d_,R_gt,t_gt]=Build_Environment(n_ele,noise,outlier_ratio,scale_gt,show_figure);

tic;

solved=0;

% for trials=1:2

scale_set=[];

break_con=1;

break_num_con=3;

% if trials==1
noise_bound=3.63;noise_bound_=5.4;
noise_bound=4.3;noise_bound_=5.2;

% else
% noise_bound=3.11;noise_bound_=4.8; 
% end

for itr_RANSAC=1:1e+9

scale_set_this=randperm(n_ele,3);

%compute controids and scale-based invariants

q_=zeros(3,1);
p_=zeros(3,1);

for i=1:3

q_(1)=q_(1)+pts_3d_(scale_set_this(i),1);
q_(2)=q_(2)+pts_3d_(scale_set_this(i),2);
q_(3)=q_(3)+pts_3d_(scale_set_this(i),3);

p_(1)=p_(1)+pts_3d(scale_set_this(i),1);
p_(2)=p_(2)+pts_3d(scale_set_this(i),2);
p_(3)=p_(3)+pts_3d(scale_set_this(i),3);

end

p_=p_/3;
q_=q_/3;

s_=zeros(3,1);s_weight=zeros(3,1);

for i=1:3
    p_tilde=norm(pts_3d(scale_set_this(i),:)'-p_);
    s_(i)=norm(pts_3d_(scale_set_this(i),:)'-q_)/p_tilde;
    s_weight(i)=p_tilde^2/0.01;
end

    mean_s=(s_(1)*s_weight(1)+s_(2)*s_weight(2)+s_(3)*s_weight(3))/(s_weight(1)+s_weight(2)+s_weight(3));

if known_scale==0
    Comp=1;
    s__=mean_s;
else
    s__=1;
    Comp=[abs(mean_s-1)<=0.06];
end


if   max(s_)-min(s_)<=0.1*s__ && Comp && ...
     abs(s_(1)-s_(2))<=noise_bound*noise/norm(pts_3d(scale_set_this(1),:)'-p_)+noise_bound*noise/norm(pts_3d(scale_set_this(2),:)'-p_) && ...
     abs(s_(1)-s_(3))<=noise_bound*noise/norm(pts_3d(scale_set_this(1),:)'-p_)+noise_bound*noise/norm(pts_3d(scale_set_this(3),:)'-p_) && ...
     abs(s_(2)-s_(3))<=noise_bound*noise/norm(pts_3d(scale_set_this(2),:)'-p_)+noise_bound*noise/norm(pts_3d(scale_set_this(3),:)'-p_) %5



if  abs(norm(pts_3d_(scale_set_this(1),:)-pts_3d_(scale_set_this(2),:))-mean_s*norm(pts_3d(scale_set_this(1),:)-pts_3d(scale_set_this(2),:)))<=2*noise*noise_bound && ...
    abs(norm(pts_3d_(scale_set_this(1),:)-pts_3d_(scale_set_this(3),:))-mean_s*norm(pts_3d(scale_set_this(1),:)-pts_3d(scale_set_this(3),:)))<=2*noise*noise_bound && ...
    abs(norm(pts_3d_(scale_set_this(2),:)-pts_3d_(scale_set_this(3),:))-mean_s*norm(pts_3d(scale_set_this(2),:)-pts_3d(scale_set_this(3),:)))<=2*noise*noise_bound %5
    
 
    v12=pts_3d(scale_set_this(2),:)-pts_3d(scale_set_this(1),:);
    X_axis=v12'/norm(v12);
    v13=pts_3d(scale_set_this(3),:)-pts_3d(scale_set_this(1),:);
    v23=cross(v12,v13);
    Y_axis=v23'/norm(v23);
    Z_axis=cross(X_axis,Y_axis);
    
    v12=pts_3d_(scale_set_this(2),:)-pts_3d_(scale_set_this(1),:);
    X_axis_=v12'/norm(v12);
    v13=pts_3d_(scale_set_this(3),:)-pts_3d_(scale_set_this(1),:);
    v23=cross(v12,v13);
    Y_axis_=v23'/norm(v23);
    Z_axis_=cross(X_axis_,Y_axis_);
    
    R_raw=[X_axis_,Y_axis_,Z_axis_]*[X_axis,Y_axis,Z_axis]';
    
    t_raw1=(pts_3d_(scale_set_this(1),:))'-mean_s*R_raw*((pts_3d(scale_set_this(1),:)))';
    t_raw2=(pts_3d_(scale_set_this(2),:))'-mean_s*R_raw*((pts_3d(scale_set_this(2),:)))';
    t_raw3=(pts_3d_(scale_set_this(3),:))'-mean_s*R_raw*((pts_3d(scale_set_this(3),:)))';
    
    check_increase=0;
    
    if  norm(t_raw1 - t_raw2) <= noise*noise_bound_ && norm(t_raw1 - t_raw3) <= noise*noise_bound_ && ...
        norm(t_raw2 - t_raw3) <= noise*noise_bound_ 

    scale_set=[scale_set;[scale_set_this,mean_s,reshape(R_raw,[9,1])',t_raw1',t_raw2',t_raw3',norm(scale_set_this)]];
    
    check_increase=1;
 
    end

    break_it=0;
    if size(scale_set,1)>=2 && check_increase==1
      H_mat=[];
      for iii=1:size(scale_set,1)-1
          
         check_scale=abs(scale_set(end,4)-scale_set(iii,4))/mean([scale_set(end,4),scale_set(iii,4)]);
         
         if check_scale~=0 &&  scale_set(end,end)~=scale_set(iii,end) && ...
            check_scale <=0.1
            
            check_rot=AngularError(reshape(scale_set(end,5:13),[3,3]),reshape(scale_set(iii,5:13),[3,3]))*180/pi;
            check_tran=max([norm(scale_set(end,14:16) - scale_set(iii,14:16)),...
             norm(scale_set(end,17:19) - scale_set(iii,17:19)), ...
             norm(scale_set(end,20:22) - scale_set(iii,20:22)),...
             norm(scale_set(end,14:16) - scale_set(iii,17:19)),...
             norm(scale_set(end,17:19) - scale_set(iii,20:22)),...
             norm(scale_set(end,14:16) - scale_set(iii,20:22))]);
         
            if check_tran <= 2*noise*5.2 && check_rot <= 10
              break_it=break_it+1;
              H_mat=[H_mat,iii];
            end
         end
      end
      H_mat=[H_mat,size(scale_set,1)];
    end
    
    
    if break_it>=break_con
        num_in=length(unique(scale_set(H_mat,1:3)));
        if num_in>=break_num_con
            
inlier_set=unique(reshape(scale_set(H_mat,1:3),1,[]));

scale_set_this=inlier_set;

q_=zeros(3,1);
p_=zeros(3,1);

len=length(inlier_set);

for i=1:len

q_(1)=q_(1)+pts_3d_(scale_set_this(i),1);
q_(2)=q_(2)+pts_3d_(scale_set_this(i),2);
q_(3)=q_(3)+pts_3d_(scale_set_this(i),3);

p_(1)=p_(1)+pts_3d(scale_set_this(i),1);
p_(2)=p_(2)+pts_3d(scale_set_this(i),2);
p_(3)=p_(3)+pts_3d(scale_set_this(i),3);

end

p_=p_/len;
q_=q_/len;

s_=zeros(len,1);s_weight=zeros(len,1);p_tilde=zeros(len,1);

for i=1:len
    p_tilde(i)=norm(pts_3d(scale_set_this(i),:)'-p_);
    s_(i)=norm(pts_3d_(scale_set_this(i),:)'-q_)/p_tilde(i);
    s_weight(i)=p_tilde(i)^2/0.01;
end

scale_opt=0;

for i=1:len
scale_opt=scale_opt+s_(i)*s_weight(i);
end

scale_opt=scale_opt/sum(s_weight);


scale_set_output=scale_set(H_mat,:);
            
n_ele_=length(inlier_set);

pts_3d_new=pts_3d(inlier_set,:);

pts_3d_new_=pts_3d_(inlier_set,:);


H=zeros(3,3);
for i=1:n_ele_
    H=H+(pts_3d_new(i,:)'-p_)*(pts_3d_new_(i,:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';
            
t_opt=q_-scale_opt*R_opt*p_;

inlier_set=zeros(1,1);
outlier_error=zeros(1,1);
inlier_error=zeros(1,1);
coun=0;
counn=0;
for i=1:n_ele
   error_this=norm(scale_opt*R_opt*pts_3d(i,:)'+t_opt-pts_3d_(i,:)');
    if error_this<=noise*1.7*3
        coun=coun+1;
        inlier_set(coun)=i;
        inlier_error(coun)=error_this;
    else counn=counn+1;
        outlier_error(counn)=error_this;
    end
end

if length(inlier_error)>=9 && sqrt(inlier_error*inlier_error'/coun)<=2.78*noise %mean(inlier_error)<=2*noise %1.73*noise*1.8 
    solved=1;
    break
end

% if trials==1
%     if break_con>=3 %|| size(scale_set,1)>=500
%         break
%     end
% else
    if break_con>=3
        break
    end
% end

break_con=break_con+1;

if break_con>=3
   break_num_con=9;
else
   break_num_con=5;
end
            
        end
    end
  
end

end

end

% if solved==1
%     break
% end
% 
% end

n_ele_=length(inlier_set);

pts_3d_new=pts_3d(inlier_set,:);

pts_3d_new_=pts_3d_(inlier_set,:);

q_=zeros(3,1);
p_=zeros(3,1);

for i=1:n_ele_

q_(1)=q_(1)+pts_3d_new_(i,1);
q_(2)=q_(2)+pts_3d_new_(i,2);
q_(3)=q_(3)+pts_3d_new_(i,3);

p_(1)=p_(1)+pts_3d_new(i,1);
p_(2)=p_(2)+pts_3d_new(i,2);
p_(3)=p_(3)+pts_3d_new(i,3);

end

p_=p_/n_ele_;
q_=q_/n_ele_;

s_=zeros(n_ele_,1);s_weight=zeros(n_ele_,1);p_tilde=zeros(n_ele_,1);

for i=1:n_ele_
    p_tilde(i)=norm(pts_3d_new(i,:)'-p_);
    s_(i)=norm(pts_3d_new_(i,:)'-q_)/p_tilde(i);
    s_weight(i)=p_tilde(i)^2/0.01;
end

scale_opt=0;

for i=1:n_ele_
scale_opt=scale_opt+s_(i)*s_weight(i);
end

scale_opt=scale_opt/sum(s_weight);

H=zeros(3,3);
for i=1:n_ele_
    H=H+(pts_3d_new(i,:)'-p_)*(pts_3d_new_(i,:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-scale_opt*R_opt*p_;

time=toc();

R_error=AngularError(R_gt,R_opt)*180/pi;
t_error=norm(t_opt - t_gt');
s_error=abs(scale_opt - scale_gt);
recall=min([100,length(inlier_set)/(n_ele*(1-outlier_ratio))*100]);


%% Display Results: 

disp(['Ground-truth Scale : ']);

scale_gt

disp(['Estimated Scale : ']);

scale_opt

disp(['Ground-truth Rotation : ']);

R_gt

disp(['Estimated Rotation : ']);

R_opt

disp(['Ground-truth Translation : ']);

t_gt'

disp(['Estimated Translation : ']);

t_opt

disp(['Scale Error : ', num2str(s_error)]);

disp(['Rotation Error (in degree): ', num2str(R_error)]);

disp(['Translation Error : ', num2str(t_error)]);

disp(['Runtime of RANSIC (in second): ', num2str(time)]);

disp(['Recall Ratio (in percentage): ', num2str(recall)]);

