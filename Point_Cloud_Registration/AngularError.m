function err = AngularError(R_pre, R_gt)

err = abs(acos((trace(R_pre'*R_gt)-1)/2));

end