function [hist_features size_hist] = hist_FE(hist_consts,...
                                                    max_hist_size,frame)


    max_hist_amp=hist_consts.max_hist_amp;
    frame_filt=frame;


    hist_from_frames_X_filt=sum(frame_filt,1);
    hist_from_frames_Y_filt=sum(frame_filt,2);

    size_hist=[0 0];
    hist_features = zeros(2*max_hist_size,1);
    if (min(length(hist_from_frames_X_filt),...
            length(hist_from_frames_Y_filt))>2)
        

        len_histX=length(hist_from_frames_X_filt);
        len_histY=length(hist_from_frames_Y_filt);
        len_hist=max(len_histX,len_histY);
        size_hist=[len_histX len_histY];

        
        if (len_hist>3)
            
            xhist_scaled = get_first_nhist(hist_from_frames_X_filt,...
                            max_hist_size);
            yhist_scaled = get_first_nhist(hist_from_frames_Y_filt,...
                            max_hist_size);

            factor_amp = 1;
            feat_id=1:2*max_hist_size;
            hist_features(feat_id(1:end/2))=xhist_scaled*factor_amp;
            hist_features(feat_id(1+end/2:end))=yhist_scaled*factor_amp;      
            
        end
        
    end

end