function [pt, N] = innerposition_exact(filt_bin_imagepart,x1,x2,y1,y2,ts)
threshold = 1;
num_region = 8;
x_hist = sum(filt_bin_imagepart,1);
y_hist = sum(filt_bin_imagepart,2);
[xhist_pos_region,num_xposregion]=gethist_inner(num_region,x_hist,x1,x2);
[yhist_pos_region,num_yposregion]=gethist_inner(num_region,y_hist,y1,y2);

N = 0;
pt = [];
if num_xposregion>0 && num_yposregion>0
   for j = 1:num_xposregion
       for k = 1:num_yposregion
           x_temp1 = xhist_pos_region(2*j-1);
           x_temp2 = xhist_pos_region(2*j);
           y_temp1 = yhist_pos_region(2*k-1);
           y_temp2 = yhist_pos_region(2*k);
           
           cnt = sum(sum(filt_bin_imagepart(y_temp1:y_temp2,x_temp1:x_temp2)));
           
           if cnt>threshold
              %if (y_temp2-y_temp1+1)*(x_temp2-x_temp1+1)<160000
                  pt = annotation_point(x1+x_temp1-1, y1+y_temp1-1, x_temp2-x_temp1+1, y_temp2-y_temp1+1, ts);
                  N = N+1;
                  if (N==num_region)
                      break;
                  end
              %end
           end
       end
       if (N==num_region)
           break;
       end
   end
end