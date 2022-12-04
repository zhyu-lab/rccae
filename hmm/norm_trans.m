function T = norm_trans(T,clamp_thres)
%this function is based on MK_STOCHASTIC that ensures the sum over every row
% of the matrix is 1. Also clamp_thres is used to make sure that
%the diagonal elements are no less than it.

for i = 1:size(T,1)
    temp = T(i,:);
    tmp = sum(temp);
    if T(i,i)<clamp_thres*tmp 
        temp(i) = 0;
        T(i,:) = (1-clamp_thres).*temp/sum(temp);
        T(i,i) = clamp_thres;
    else
        T(i,:) = temp/(tmp+eps);%to avoid tmp = 0
    end
end



