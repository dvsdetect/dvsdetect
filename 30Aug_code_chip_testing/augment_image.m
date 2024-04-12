function [aug_hist_feat, n_img] = augment_image(hist_consts,...
	                              image,fixed_size,...
	                              features)
    [h,w] = size(image);
    aug_hist_feat = [];
    rng(0,'twister');
    function aug_hist = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist)
		[histf, histf_size] = hist_FE(hist_consts,...
	                                    fixed_size,inp_img);
	    aug_hist = [aug_hist,histf ];
	end
    if (w>fixed_size) && (h>fixed_size)

        if strcmp(features,"hist_FE")
        	inp_img = image(1:fixed_size,(w-fixed_size+1):end);
            
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            inp_img = image((h - fixed_size+1):end,1:fixed_size);
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            inp_img = image((h - fixed_size+1):end,(w-fixed_size+1):end);
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            n_img = 3;
        end
        
    elseif (w>fixed_size) && (h<=fixed_size)
        if strcmp(features,"hist_FE")
        	inp_img = image(:,(w-fixed_size+1):end);
            
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            n_img = 1;
        end
    elseif (h>fixed_size) && (w<=fixed_size)
        if strcmp(features,"hist_FE")
        
            inp_img = image((h - fixed_size+1):end,:);
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            n_img = 1;
        end
    else
    	if strcmp(features,"hist_FE")
            a = -15;
            b = 15;
            r = (b-a).*rand(2,1) + a;
            inp_img = imrotate(image,r(1),'bilinear','crop');
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            inp_img = imrotate(image,r(2),'bilinear','crop');
            aug_hist_feat = add_feat(hist_consts,fixed_size,...
							inp_img,aug_hist_feat);
            n_img = 2;
        end
    end
    	
end

