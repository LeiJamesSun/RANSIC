%% Solver 'RANSIC' for Rotation Search

% This is the source code of paper: 'RANSIC: Fast and Highly Robust 
% Estimation for Rotation Search and Point Cloud Registration using 
% Invariant Compatibility'

% Copyright by Lei Sun (leisunjames@126.com)

clc;
clear all;
close all;

%% Parameter Setup by User(s): 

% set correspondence number N 
n_ele=100; % should be set to 100 or 500 or 1000

% set outlier ratio in percentage
outlier_ratio=0.90; % should be set within 0-0.99 (for successful results, 0-0.95 with N=100, 0-0.98 with N=500, 0-0.99 with N=1000)



%% RANSIC Begins: (No Need for Adjustment)

disp(['RANSIC for Rotation Search :']);

disp(['N=',num2str(n_ele)]);

%set noise in pixels
noise=0.01;

disp(['Noise=',num2str(noise)]);

% build simulation environment

% Input:
% pts_3d is a Nx3 matrix
% pts_3d_ is a Nx3 matrix
% R_gt is a 3x3 SO(3) matrix ('gt' means ground-truth)

[pts_3d,pts_3d_,R_gt]=Build_Environment(n_ele,noise,outlier_ratio);

% set bounds or conditions

if n_ele==100
    bound1=0.012;
    min_break_num_con1=4;
    min_break_num_con2=3;
    stop_con=3.2;
    if outlier_ratio<=0.95
        min_in_size=5;
        min_ang=5;
    else
        min_in_size=4;
        min_ang=6;
    end
elseif n_ele==500
    bound1=0.008;
    min_ang=5;
    stop_con=2.7;
    if outlier_ratio<=0.98
        min_in_size=10;
        min_break_num_con1=8;
        min_break_num_con2=5;
    else
        min_in_size=5;
        min_break_num_con1=4;
        min_break_num_con2=3;
    end
elseif n_ele==1000
    bound1=0.008;
    min_ang=5;
    min_in_size=9;
    min_break_num_con1=8;
    min_break_num_con2=5;
    stop_con=2.7;
end


tic;

scale_set=[];

break_con=1;

break_num_con=3;

for itr_RANSAC=1:1e+8

scale_set_this=randperm(n_ele,2);


length_a=norm(pts_3d(scale_set_this(1),:)-pts_3d(scale_set_this(2),:));

length_b=norm(pts_3d_(scale_set_this(1),:)-pts_3d_(scale_set_this(2),:));

check_increase=0;

if abs((length_a-length_b))<=bound1
    
    v12=pts_3d(scale_set_this(1),:);
    X_axis=v12';
    v13=pts_3d(scale_set_this(2),:);
    v23=cross(v12,v13);
    Y_axis=v23'/norm(v23);
    Z_axis=cross(X_axis,Y_axis);
    
    v12=pts_3d_(scale_set_this(1),:);
    X_axis_=v12';
    v13=pts_3d_(scale_set_this(2),:);
    v23=cross(v12,v13);
    Y_axis_=v23'/norm(v23);
    Z_axis_=cross(X_axis_,Y_axis_);
    
    R_raw=[X_axis_,Y_axis_,Z_axis_]*[X_axis,Y_axis,Z_axis]';
    
    scale_set=[scale_set;[scale_set_this,reshape(R_raw,[9,1])',norm(scale_set_this)]];
    
    check_increase=1;
    
end

    break_it=0;
    if size(scale_set,1)>=2 && check_increase==1
      H_mat=[];
      for iii=1:size(scale_set,1)-1
          
            
            check_rot=AngularError(reshape(scale_set(end,3:11),[3,3]),reshape(scale_set(iii,3:11),[3,3]))*180/pi;
         
            if check_rot<= min_ang
              break_it=break_it+1;
              H_mat=[H_mat,iii];
            end
      end
      H_mat=[H_mat,size(scale_set,1)];
    end
    
    
    if break_it>=break_con
        num_in=length(unique(scale_set(H_mat,1:2)));
        if num_in>=break_num_con
            
inlier_set=unique(reshape(scale_set(H_mat,1:2),1,[]));

H=zeros(3,3);
for i=1:length(inlier_set)
    H=H+(pts_3d(inlier_set(i),:)')*pts_3d_(inlier_set(i),:);
end

[U,~,V]=svd(H);

R_opt=V*U';

inlier_set=zeros(1,1);coun=0;counn=0;res_in=zeros(1,1);res_out=zeros(1,1);

for i=1:n_ele
   
    res_=norm(R_opt*pts_3d(i,:)'-pts_3d_(i,:)');
    
    if res_<=noise*1.73*3
        coun=coun+1;
        inlier_set(coun)=i;
        res_in(coun)=res_;
    else counn=counn+1;
        res_out(counn)=res_;
    end
    
end
        
if length(res_in)>=min_in_size && sqrt(res_in*res_in'/coun)<=stop_con*noise %mean(res_in)<=noise*1.73*1.5 
    break
end

if break_con>=3
    break
end

break_con=break_con+1;

if break_con>=3
   break_num_con=min_break_num_con2;
else
   break_num_con=min_break_num_con1;
end
         
            
        end
    end
  
end


n_ele_=length(inlier_set);

pts_3d_new=pts_3d(inlier_set,:);

pts_3d_new_=pts_3d_(inlier_set,:);

H=zeros(3,3);
for i=1:n_ele_
    H=H+(pts_3d_new(i,:)')*pts_3d_new_(i,:);
end

[U,~,V]=svd(H);

R_opt=V*U';

R_error=AngularError(R_gt,R_opt)*180/pi;

time=toc();

Recall=min([100,length(inlier_set)/(n_ele*(1-outlier_ratio))*100]);


%% Display Results: 

disp(['Ground-truth Rotation : ']);

R_gt

disp(['Estimated Rotation : ']);

R_opt

disp(['Rotation Error (in degree): ', num2str(R_error)]);

disp(['Runtime of RANSIC (in second): ', num2str(time)]);

disp(['Recall Ratio (in percentage): ', num2str(Recall)]);
