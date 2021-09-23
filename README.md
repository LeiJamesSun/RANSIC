# RANSIC
'RANSIC: Fast and Highly Robust Estimation for Rotation Search and Point Cloud Registration using Invariant Compatibility'

Copyright by Lei Sun (leisunjames@126.com)

This package is the source code of our solver 'RANSIC' in Matlab, only for academic use.
 
%%%%%%%%%%

(A) Standard experiments:

To implement RANSIC for rotation search, please:
1. open 'RANSIC_for_RS.m' in the file 'Rotation_Search';
2. adjust the parameters (correspondence number, outlier ratio) as suggested in the comments;
3. run the script 'RANSIC_for_RS.m' for results.

To implement RANSIC for point cloud registration, please:
1. open 'RANSIC_for_PCR.m' in the file 'Point_Cloud_Registration';
2. adjust the parameters (scale, outlier ratio) as suggested in the comments;
3. run the script 'RANSIC_for_PCR.m' for results.

Note that the details of the parameter setup are specified in the comments.


(B) Real applications:

To test RANSIC for image stitching, please:
1. open 'Demo_Img_Stitching.m' in the file 'Rotation_Search';
2.choose whether to show the stitching result according to the comments (depending on the CPU and RAM);
3.run the script for both qualitative and quantitative results.

To test RANSIC for object localization, please:
1. open 'Demo_Obj_Localization.m' in the file 'Point_Cloud_Registration';
2.run the script for both qualitative and quantitative results.

%%%%%%%%%%

Please have a good time implementing our solver!


