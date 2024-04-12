function [pt, N] = sepposition_exact5slens(filt_bin_image,filt_bin_imagepart,x1,x2,y1,y2,ts)
a = filt_bin_imagepart;
p1 = find(all(a==0,2));
p2 = find(all(a==0,1));

N = 0;
pt = [];
threshold = 5;
num_region = 8;
gap = 15;

if length(p1)>1 && length(p2)>1
   x_hist = sum(filt_bin_imagepart(p1(1):p1(end),p2(1):p2(end)),1);
   y_hist = sum(filt_bin_imagepart(p1(1):p1(end),p2(1):p2(end)),2);
   [xhist_pos_region,num_xposregion]=gethist_inner(num_region,x_hist,p2(1),p2(end));
   [yhist_pos_region,num_yposregion]=gethist_inner(num_region,y_hist,p1(1),p1(end));

   if num_xposregion>0 && num_yposregion>0
      for j = 1:num_xposregion
          for k = 1:num_yposregion
              x_temp1 = xhist_pos_region(2*j-1);
              x_temp2 = xhist_pos_region(2*j);
              y_temp1 = yhist_pos_region(2*k-1);
              y_temp2 = yhist_pos_region(2*k);
           
               cnt = sum(sum(filt_bin_imagepart(y_temp1:y_temp2,x_temp1:x_temp2)));
               cnt2 = sum(sum(filt_bin_image(max(1,y1+y_temp1-gap):min(260,y1+y_temp2+gap),max(1,x1+x_temp1-gap):min(346,x1+x_temp2+gap))));
               cnt3 = sum(sum(filt_bin_image(1:min(260,y1+y_temp1-2),max(1,x1+x_temp1):min(346,x1+x_temp2))));
           
              if cnt2-cnt<=3 && cnt3<=3 && cnt>threshold && (x_temp2-x_temp1+1)*(y_temp2-y_temp1+1)> 5*5 && cnt/((x_temp2-x_temp1+1)*(y_temp2-y_temp1+1))>0.3 && cnt/((x_temp2-x_temp1+1)*(y_temp2-y_temp1+1))<0.95
                 pt = annotation_point(x1+x_temp1-1, y1+y_temp1-1, x_temp2-x_temp1+1, y_temp2-y_temp1+1, ts);
                 N = N+1;
                 if (N==num_region)
                    break;
                 end
              end
          end
          if (N==num_region)
              break;
          end
      end
   end
end