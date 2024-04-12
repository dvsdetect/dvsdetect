function [xhist_pos_region,num_xposregion]=gethist(num_region,a_maxXSize,scale,threshold,xhist_pos_s4)
num_xposregion=0;
start=0;
xhist_pos_region=zeros(2*num_region,1);
for j=1:floor(a_maxXSize/scale) %find regions in scaled X histogram,a_maxXSize/4=length(xhist_s4)
    if xhist_pos_s4(j)>threshold 
        if start==0
            start=1;
            num_xposregion=num_xposregion+1;
            if num_xposregion>num_region
                num_xposregion=num_region;
                return;
            end
            xhist_pos_region(2*num_xposregion-1)=scale*j-scale+1;
        end
    else
       if start==1
           xhist_pos_region(2*num_xposregion)=scale*(j-1);
           start=0;
       end
    end
end
if start==1
    xhist_pos_region(2*num_xposregion)=a_maxXSize;
    start=0;
end