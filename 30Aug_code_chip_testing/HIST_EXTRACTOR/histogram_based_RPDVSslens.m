function [boxes, N] = histogram_based_RPDVSslens(consts,...
                        filt_bin_image,ts)
    num_region = consts.RP_consts.n_objs;        
    threshold = consts.RP_consts.th;
    scale = 2;
    a_maxXSize = consts.maxX;
    a_maxYSize = consts.maxY;
    xhist       = sum(filt_bin_image,1);
    yhist       = sum(filt_bin_image,2);
    
    yhist_s1    = sum(reshape(yhist,scale,floor(a_maxYSize/scale)));
    xhist_s2    = sum(reshape(xhist,scale,floor(a_maxXSize/scale)));
    [xhist_region,num_xregion]=gethist(num_region,a_maxXSize,...
                              scale,threshold,xhist_s2);
    [yhist_region,num_yregion]=gethist(num_region,a_maxYSize,...
                               scale,threshold,yhist_s1);
    

    boxes = [];
    N = 0;
    Num_inner = 0;
    Num_inner1 = 0;
    Num_inner2 = 0;
    
    pt_temp = [];
    
    if num_xregion >0 && num_yregion >0
        for j=1:num_xregion
            for k=1:num_yregion
                x1=xhist_region(2*j-1);
                x2=xhist_region(2*j);
                y1=yhist_region(2*k-1);
                y2=yhist_region(2*k);


                cnt=sum(sum(filt_bin_image(y1:y2,x1:x2)));

                if cnt>threshold

                     if (((y2-y1)*(x2-x1)<a_maxYSize*a_maxXSize/2) && ...
                         (y2-y1+1>3) && (x2-x1+1>3)) 
                     if (((y2-y1+1)*(x2-x1+1)<80*80) && ((y2-y1+1)*(x2-x1+1)>3*3) && ...
                        (y2-y1+1>3) && (x2-x1+1>3))                       
                        if y1-3<1
                           t_y1=1;
                        else
                           t_y1=y1-3;
                        end
                        if y2+3>=a_maxYSize
                           t_y2=a_maxYSize;
                        else
                           t_y2=y2+3;
                        end
                        if x1-3<=1
                           t_x1=1;
                        else
                           t_x1=x1-3;
                        end
                        if x2+3>=a_maxXSize
                           t_x2=a_maxXSize;
                        else
                           t_x2=x2+3;
                        end
                        [pt, Num_inner] = sepposition_exact5slens(filt_bin_image,filt_bin_image(t_y1:t_y2,t_x1:t_x2),t_x1,t_x2,t_y1,t_y2,ts);                        
                        boxes = [boxes, pt];
                        N = N + Num_inner;
                     end
                     if (((y2-y1+1)*(x2-x1+1)>=80*80) && ...
                          (y2-y1+1>40) && (x2-x1+1>40)) 
                        [pt1, Num_inner1] = simplifedccl_exactslens(filt_bin_image(y1:y2,x1:x2),x1,x2,y1,y2,ts);
                        boxes = [boxes, pt1];
                        N = N + Num_inner1;
                     end
                        if (N==num_region)
                            break;
                        end
                 end
             end
         end
            if (N==num_region)
                return;
            end
        end
end