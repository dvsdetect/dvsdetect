function hist_out = scale_hist(hist_orig,hist_target,factor)

max_size=length(hist_target);

for i=1:max_size
    n=i*factor;
    index=ceil(n);
    if index>length(hist_orig)
        break;
    end
    current=n-floor(n);
    hist_target(i)=hist_target(i)+current*hist_orig(index);
    remaining=factor-current;
    while remaining>1
        current=1;
        hist_target(i)=hist_target(i)+current*hist_orig(index);
        index=index-1;
        remaining=remaining-current;
    end
    hist_target(i)=hist_target(i)+remaining*hist_orig(index);

end


hist_out=hist_target;