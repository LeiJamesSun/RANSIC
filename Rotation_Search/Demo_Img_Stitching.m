%% Solver 'RANSIC' for Image Stitching

% This is the source code of paper: 'RANSIC: Fast and Highly Robust 
% Estimation for Rotation Search and Point Cloud Registration using 
% Invariant Compatibility'

% Copyright by Lei Sun (leisunjames@126.com)

clc;
clear all;
close all;

addpath include;
load('include/52_63_data.mat');

I1 = rgb2gray(imread('include/img52.jpg'));
I2 = rgb2gray(imread('include/img63.jpg'));


U1=matchedPoints1.Location;

U2=matchedPoints2.Location;

pts_3d=inv(T1)*[U1,ones(matchedPoints1.Count,1)]';

pts_3d_=inv(T2)*[U2,ones(matchedPoints2.Count,1)]';

pts_3d=pts_3d';

pts_3d_=pts_3d_';

cc=matchedPoints1.Count;

for i=1:cc

pts_3d(i,:)=double(pts_3d(i,:)/norm(pts_3d(i,:)));

pts_3d_(i,:)=double(pts_3d_(i,:)/norm(pts_3d_(i,:)));

end

pts_3d=double(pts_3d);pts_3d_=double(pts_3d_);

n_ele=cc;


%% RANSIC for Solving Rotation

noise=0.02;

for repeat=1:2

tic;

scale_set=[];

break_con=1;

break_num_con=3;


for itr_RANSAC=1:1e+8

scale_set_this=randperm(n_ele,2);


length_a=norm(pts_3d(scale_set_this(1),:)-pts_3d(scale_set_this(2),:));

length_b=norm(pts_3d_(scale_set_this(1),:)-pts_3d_(scale_set_this(2),:));

check_increase=0;

if abs((length_a-length_b))<=2*0.01 %0.008 0.1
    
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
    
    scale_set=[scale_set;[scale_set_this,Vec(R_raw)',norm(scale_set_this)]];
    
    check_increase=1;
    
end

    break_it=0;
    if size(scale_set,1)>=2 && check_increase==1
      H_mat=[];
      for iii=1:size(scale_set,1)-1
          
            
            check_rot=AngularError(reshape(scale_set(end,3:11),[3,3]),reshape(scale_set(iii,3:11),[3,3]))*180/pi;
         
            if check_rot<= 1
              break_it=break_it+1;
              H_mat=[H_mat,iii];
            end
      end
      H_mat=[H_mat,size(scale_set,1)];
    end
    
    
    if break_it>=break_con
        num_in=length(unique(scale_set(H_mat,1:2)));
        if num_in>=break_num_con
            
inlier_set=unique(Vec(scale_set(H_mat,1:2)));

H=zeros(3,3);
for i=1:length(inlier_set)
    H=H+(pts_3d(inlier_set(i),:)')*pts_3d_(inlier_set(i),:);
end

[U,~,V]=svd(H);

R_opt=V*U';

inlier_set=zeros(1,1);coun=0;counn=0;res_in=zeros(1,1);res_out=zeros(1,1);

for i=1:n_ele
   
    res_=norm(R_opt*pts_3d(i,:)'-pts_3d_(i,:)');
    
    if res_<=noise*1.73*4
        coun=coun+1;
        inlier_set(coun)=i;
        res_in(coun)=res_;
    else counn=counn+1;
        res_out(counn)=res_;
    end
    
end
        
if length(res_in)>=round((1-0.95)*n_ele) && mean(res_in)<=noise*3.2 %R_error(itr)<=5
    break
end

if break_con>=3
    break
end

break_con=break_con+1;

if break_con>=3
   break_num_con=5; %8 9(0.1) 4
else
   break_num_con=4; %5 6(0.1) 3
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

end


time=toc();

recall_RANSIC=length(inlier_set);





% Read in images

filename1='include/img52.jpg';
filename2='include/img63.jpg';

scale = 1;
I1 = imresize(imread(filename1),scale);
I2 = imresize(imread(filename2),scale);


% Image size
[h,w,~] = size(I1);

p1=U1';p2=U2';


% Infinite homography
H=T1*R_opt*inv(T1);

th=0.1;
inlrsIdx = anginlrs(pts_3d*R_opt',pts_3d_,th);



%% Plot Input Raw 2D Correspondences

figure, imshow([I1,I2])
hold on
for i=1:size(p1,2)
    p1i = p1(:,i);
    p2i = p2(:,i);
    plot([p1i(1), p2i(1)+w],[p1i(2), p2i(2) ],'-c','linewidth',1.2)
end



%% Plot Inliers and Outliers Obtained by RANSIC

figure, imshow([I1,I2])
hold on

for i=1:size(p1,2)
    p1i = p1(:,i);
    p2i = p2(:,i);
    plot([p1i(1), p2i(1)+w],[p1i(2), p2i(2) ],'-g','linewidth',1)
    if ismember(i,inlier_set)
        plot([p1i(1), p2i(1)+w],[p1i(2), p2i(2) ],'-g','linewidth',1)
    else
        plot([p1i(1), p2i(1)+w],[p1i(2), p2i(2) ],'-r','linewidth',1)
    end
end



%% Plot Blended images (optional, depending on the CPU and RAM, 32GB of RAM is recommended)

% To show the image stitching results, please uncomment the following codes: (better implemented with 32GB or more RAM)

% figure;
% I1p = paddedimage(I1,1);
% imgwarped  = image_warping(I1p,I2,H);
% imgblended = image_blending(I1p,imgwarped);
% imshow(imgblended);



%% Display Results: 

disp(['Runtime of RANSIC (in second): ', num2str(time)]);

disp(['Inliers Found by RANSIC: ', num2str(recall_RANSIC)]);

disp('The first image shows the raw correspondences matched by SURF;');

disp('The second image shows the inlier correspondences found by RANSIC;');

disp('The third image shows stitched images. (This is optional, depending on the CPU and RAM.)');

