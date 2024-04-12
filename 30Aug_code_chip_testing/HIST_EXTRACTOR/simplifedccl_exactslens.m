function [pt, N] = simplifedccl_exactslens(filt_bin_imagepart,x1,x2,y1,y2,ts)
gap = 15;
b = filt_bin_imagepart;
N = 0;
pt = [];
pt_temp = [];
num_region = 8;

total_pixel = sum(sum(b));

j = 1;
record_p1 = ones(1,size(b,2))*size(b,1);
record_p2 = ones(1,size(b,2))*size(b,1);
for i = 1:1:size(b,2)
    temp = b(:,i);
    [row1,col1] = find(temp == 1);
    if isempty(row1) == 0
       record_p1(1,i) = min(row1);
       record_p2(1,i) = i;
       j = j+1;
    end
end

j1 = 1;
startpoint = [];
pointlength = [];
buffer = [];
length = 0;
for k = 2:1:size(b,2)
       if record_p1(1,k)-record_p1(1,k-1)<(-1*gap)
          startpoint(1,j1) = k;
          j1 = j1+1;
       end
end

for ii = 1:1:size(startpoint,2)
    if ii<size(startpoint,2)
       if startpoint(1,ii+1)-startpoint(1,ii)<gap
          pointlength(ii) = 1;
       end
       if startpoint(1,ii+1)-startpoint(1,ii)>gap
          compute1 = record_p1(1,startpoint(ii):startpoint(ii+1)-2);
          compute2 = record_p1(1,startpoint(ii)+1:startpoint(ii+1)-1);
          compute21 = compute2-compute1;
          temp_compute = [0 compute21];
          [rw,cl] = find(temp_compute>gap);
          if isempty(cl)==0
             pointlength(ii) = min(cl);
          else
             pointlength(ii) = 1;
          end
       end
    end
    if ii==size(startpoint,2)
       endpos = size(b,2);
       if endpos-startpoint(1,ii)<gap
          pointlength(ii) = 1;
       end
       if endpos-startpoint(1,ii)>gap
          compute1 = record_p1(1,startpoint(ii):endpos-1);
          compute2 = record_p1(1,startpoint(ii)+1:endpos);
          compute21 = compute2-compute1;
          temp_compute = [0 compute21];
          [rw,cl] = find(temp_compute>gap);
          if isempty(cl)==0
             pointlength(ii) = min(cl);
          else
             pointlength(ii) = 1;
          end
       end
    end
end

for count = 1:1:size(pointlength,2)
    if pointlength(1,count)>5    
       col_start = startpoint(1,count);
       col_end   = startpoint(1,count)+pointlength(1,count)-1;
       temp_col1 = sum(b(:,col_start:col_end),2);
       [row2,col2] = find(temp_col1>0);              
       row_start = min(row2);
       temp_col2 = sum(b(row_start:end,col_start:col_end),2);
       [row3,col3] = find(temp_col2==0);
       row_end = row_start+min(row3)-1;
       num_pixel = sum(sum(b(row_start:row_end,col_start:col_end)));
       temp_width = col_end-col_start+1;
       temp_height = row_end-row_start+1;
       if num_pixel>5 && num_pixel<0.1*total_pixel && temp_width*temp_height>3*3   
          if num_pixel/(temp_width*temp_height)>0.4 && num_pixel/(temp_width*temp_height)<0.95
             if temp_height/temp_width<1 && temp_height/temp_width>1/4
                pt_temp = annotation_point(x1+col_start-1, y1+row_start-1, col_end-col_start+1, row_end-row_start+1, ts);
                pt = [pt pt_temp];
                N = N+1;
                if (N==num_region)
                   break;
                end
             end
          end
       end
    end
end