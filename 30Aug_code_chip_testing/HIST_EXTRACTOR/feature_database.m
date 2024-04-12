%%% GT box stores the centroid of the box rather than the corner.

classdef feature_database < handle
    properties 
        consts
        meta_data = meta_data_class;
        hist_features = []
        img_features
        current_track_ID = 0
        max_track_ID_found = 0
        files_explored = 0
        all_examples = 0
        examples_by_class
        tracks_by_class
    end
    methods(Access =public)
        function obj = feature_database(consts)
            if nargin > 0
                obj.consts = consts;
                obj.meta_data = meta_data_class;
                obj.hist_features = [];
                obj.examples_by_class = zeros(1,consts.no_classes);
                obj.tracks_by_class = zeros(1,consts.no_classes);
            end
        end
        function load_file(obj,mat_file)
            obj.hist_features = [obj.hist_features,mat_file(7:end,:)];
            obj.meta_data.load_meta(mat_file(1:6,:));
            for jj = 1:obj.consts.no_classes
                ex_for_class = sum(mat_file(4,:)==jj);
                obj.examples_by_class(jj) = obj.examples_by_class(jj) ...
                                      + ex_for_class;
                tracks_for_class = mat_file(1,mat_file(4,:)==jj);
                n_uniq_tracks = length(unique(tracks_for_class));
                obj.tracks_by_class(jj) = obj.tracks_by_class(jj)...
                                         + n_uniq_tracks;
                obj.all_examples = obj.all_examples + ex_for_class;
            end
        end
        %% save_database: function description
        function save_filename = save_database(obj,filename)
            if ~exist(obj.consts.database_save_path,...
                    'dir')
                mkdir(obj.consts.database_save_path);
            end
            formatOut = 'yyyymmdd';
            dt = datestr(now,formatOut);
            save_filename = [obj.consts.database_save_path,'/',...
                            'hist_features_',filename,'_',dt,'.mat'];
            arr_to_save = [obj.meta_data.make_mat();
                            obj.hist_features];
            save(save_filename,'arr_to_save');
            obj.hist_features = [];
            obj.meta_data = meta_data_class;
        end
        
        function remove_class(obj, classID)
            if (classID >obj.consts.no_classes ) || ...
                    (classID <1)
                error(['Class not found: ',num2str(classID)]);
            else
                class_inds = obj.meta_data.classIDs==classID;
                obj.meta_data = RemoveNulls(obj.meta_data,...
                                    class_inds);
                obj.all_examples = obj.all_examples - ...
                                    obj.examples_by_class(classID);
                obj.examples_by_class(classID) = [];
                obj.tracks_by_class(classID) = [];
                obj.hist_features(:,class_inds) = [];
                obj.consts.no_classes = obj.consts.no_classes - 1;
                if (classID ~= obj.consts.no_classes)
                    t_ind = find(obj.meta_data.classIDs>classID);
                    obj.meta_data.classIDs(t_ind) = ...
                        obj.meta_data.classIDs(t_ind) - 1;
                end
            end
        end
        
        function remove_augmented_dat(obj)
            inds = logical(obj.meta_data.augmented);
            obj.hist_features(:,inds) = [];
            obj.meta_data = RemoveNulls(obj.meta_data,...
                                    inds);
            obj.all_examples = 0;
            for jj = 1:obj.consts.no_classes
                ex_for_class = sum(obj.meta_data.classIDs==jj);
                obj.examples_by_class(jj) = ex_for_class;
                
                obj.all_examples = obj.all_examples + ex_for_class;
            end
        end
        function merge_classes(obj,classIDs)
            % format: class 1 & class 2. class 1 < class 2
            % class 2 to be merged with class 1
            if (classIDs(1) >obj.consts.no_classes ) || ...
                    (classIDs(1) <1) || (classIDs(2) >...
                    obj.consts.no_classes ) || ...
                    (classIDs(2) <1) || (classIDs(1) > classIDs(2))
                error(['Class not found: ',num2str(classIDs)]);
            else
                class_inds = find(obj.meta_data.classIDs==classIDs(2));
                obj.meta_data.classIDs(class_inds) = classIDs(1);
                obj.examples_by_class(classIDs(1)) = ...
                    obj.examples_by_class(classIDs(1)) + ...
                    obj.examples_by_class(classIDs(2));
                obj.examples_by_class(classIDs(2)) = [];
                
                obj.tracks_by_class(classIDs(1)) = ...
                    obj.tracks_by_class(classIDs(1)) + ...
                    obj.tracks_by_class(classIDs(2));
                obj.tracks_by_class(classIDs(2)) = [];
                
                obj.consts.no_classes = obj.consts.no_classes - 1;
                
                if (classIDs(2) ~= obj.consts.no_classes)
                    t_ind = find(obj.meta_data.classIDs>classIDs(2));
                    obj.meta_data.classIDs(t_ind) = ...
                        obj.meta_data.classIDs(t_ind) - 1;
                end
            end
        end
        

        
        function balance_DB(obj,examples_reqd)
            rng(0,'twister');
             
            for jj = 1:obj.consts.no_classes
                if obj.examples_by_class(jj) > examples_reqd
                    list = randperm(obj.examples_by_class(jj));
                    inds_to_remove = list((examples_reqd+1):end);
                    class_inds = find(obj.meta_data.classIDs==jj);
                    track_IDs = obj.meta_data.track_ids(...
                                class_inds(list(1:examples_reqd)));
                    inds_to_remove = class_inds(inds_to_remove);
                    templist = zeros(length(obj.meta_data.classIDs),1);
                    templist(inds_to_remove) = 1;
                    obj.hist_features(:,inds_to_remove) = [];
                    obj.meta_data = RemoveNulls(obj.meta_data,...
                                    templist);
                    obj.examples_by_class(jj) = examples_reqd;
                    
                    
                    
                    n_uniq_tracks = length(unique(track_IDs));
                    obj.tracks_by_class(jj) = n_uniq_tracks;
                else
                    to_add  = examples_reqd - obj.examples_by_class(jj);
                    cycles = floor(to_add/obj.examples_by_class(jj));
                    class_inds = find(obj.meta_data.classIDs==jj);
                    added = 0;
                    for kk = 1:cycles
                        list = randperm(obj.examples_by_class(jj));
                        inds = class_inds(list);
                        obj.hist_features = [obj.hist_features , ...
                                             obj.hist_features(:,inds)];
                        obj.meta_data.copy_add(inds);
                        added  = added + obj.examples_by_class(jj);
                    end
                    to_add = examples_reqd - added- obj.examples_by_class(jj);
                    if (to_add >0)
                        list = randperm(obj.examples_by_class(jj));
                        inds = class_inds(list(1:to_add));
                        obj.hist_features = [obj.hist_features , ...
                                             obj.hist_features(:,inds)];
                        obj.meta_data.copy_add(inds);
                    end
                    obj.examples_by_class(jj) = examples_reqd;
                end
            end
        end
        
        function obj = save_RP_annotations(obj,TD,annots,...
                                        RP_type,ID)

            max_y = obj.consts.maxY;
            max_x = obj.consts.maxX;
            max_trk_i = 0;
            overlap_ratio = obj.consts.RP_consts.IoU;
            gt_classes = [annots.tracks(:).class];
            total_tracks = length(annots.tracks);
            tstart=TD.ts(1);
            
            frame_66ms = single(zeros(max_y, max_x));
            annot_writer = write_annotation(obj.consts,ID);
            ctr = 1;
            for jj = 1:length(TD.ts)
                t = TD.ts(jj);
                x = TD.x(jj);
                y = TD.y(jj);
                frame_66ms(y,x) = 1;
                
                if (t-tstart)>obj.consts.window_period
                    
                    filt_frame=...
                                my_medfilt2(frame_66ms,...
                                obj.consts.stride); 
                    if strcmp(RP_type, "histogram_based_RP")
                        [tracker_boxes, N_boxes] = ...
                                histogram_based_RP(obj.consts,...
                                               filt_frame,t);
                    elseif strcmp(RP_type, "CCA_based_RP")
                        
                        [tracker_boxes, N_boxes] = ...
                                Single_Pass_CCL(obj.consts,...
                                               filt_frame,t);
                    end

                    if (N_boxes > 0)
                        [gt_boxes, gt_trackIDs] = ...
                                annots.get_all_bounding_boxes(t);
                        cpp_classes = [];
                        cpp_trackIDs = [];
%                         disp("Entered here");
                       
                        if (length(gt_boxes)>0)
                            result = zeros(length(gt_boxes),...
                                 N_boxes); %create placeholder
                            for gt_index = 1:length(gt_boxes)
                                for cpp_index = 1:N_boxes

                                    [score1, score2] = ...
                                        compute_overlap_score2(...
                                        gt_boxes(gt_index),...
                                         tracker_boxes(cpp_index)); 
                                    % score 1 is intersection/union
                                    result(gt_index, cpp_index) = score1; 
                
                                end
                            end
                            [eval_overlaps, gt_indices] = max(result, [], 1);
                            valid_overlaps = eval_overlaps>overlap_ratio;
                            eval_overlaps = eval_overlaps(valid_overlaps); 

                            assigned_class = gt_classes(gt_trackIDs);
                            cpp_classes = assigned_class(gt_indices);
                            cpp_classes(~valid_overlaps) = 7;
                            cpp_trackIDs = gt_trackIDs(gt_indices);
                            cpp_trackIDs  =cpp_trackIDs + ...
                                            obj.max_track_ID_found;
                            cpp_trackIDs(~valid_overlaps) = total_tracks+ ctr;
                            max_trk_i = max(max_trk_i,...
                                        max(cpp_trackIDs));
                        end

                        if isempty(cpp_classes)
                            cpp_classes = 7*ones(1,N_boxes);
                            cpp_trackIDs = ones(1,N_boxes)*(total_tracks+ctr);
                        end
                        annot_writer.add_annotations(tracker_boxes,...
                            N_boxes,cpp_classes,cpp_trackIDs)
                    end
                    frame_66ms = single(zeros(max_y, max_x));
                    tstart=TD.ts(jj);  
                    ctr = ctr + 1;
                    clear tracker_boxes cpp_classes cpp_trackIDs ...
                        N_boxes 
                end
                
            end 

            annot_writer.close();   
        end
        
        function save_comparison_video(obj,TD,annots)
            
            max_y = obj.consts.maxY;
            max_x = obj.consts.maxX;
            
            
            tstart=TD.ts(1);
            
            frame_66ms = single(zeros(max_y, max_x));
            writerObj = VideoWriter('file_video.avi', 'Uncompressed AVI');
             writerObj.FrameRate = 24;
            open(writerObj);
            for jj = 1:length(TD.ts)
                t = TD.ts(jj);
                x = TD.x(jj);
                y = TD.y(jj);
                frame_66ms(y,x) = 1;
                
                if (t-tstart)>obj.consts.window_period
                    
                    filt_frame=...
                                my_medfilt2(frame_66ms,...
                                obj.consts.stride); 
                    
                    [tracker_boxes, N_boxes] = ...
                            histogram_based_RP(obj.consts,...
                                           filt_frame,t);

                    [tracker_boxes2, N_boxes2] = ...
                            Single_Pass_CCL(obj.consts,...
                                           filt_frame,t);

                    [gt_boxes, gt_trackIDs] = ...
                            annots.get_all_bounding_boxes(t);


                    if (length(gt_boxes)>0)

                        RGB = show_boxes(filt_frame,gt_boxes,...
                                    tracker_boxes,...
                                    tracker_boxes2);
                        writeVideo(writerObj,RGB);

 
 
                    end
                    
                    frame_66ms = single(zeros(max_y, max_x));
                    tstart=TD.ts(jj);  
                end
                
            end
            close(writerObj);
        end

  
          
          function [track,current] = explore_video_using_RPCCLcombineDVSslens(obj,TD,RP_type)            
            
            max_y = obj.consts.maxY;  
            max_x = obj.consts.maxX;  
            max_trk_i = 0;
            overlap_ratio = obj.consts.RP_consts.IoU;
            
            tstart=TD.ts(1);
            
            frame_temp = single(zeros(346,260));
            
            frame_ctr = 0;
            
            defnum_region = 8;
            
            
            rng(2); 

            a_maxXSize=346;                                                            
            a_maxYSize=260;                                                             

            classifier_type=1;                                                          
            numframe_per_feature=1;
            num_region=8;
            num_tracker=8;
            matched=zeros(num_region,num_tracker);                                      

            Sim_th=6;                                                                   
            history_th=2;                                                               
            classify_th=2;
            overlap_th=0.1;                                                             
            distance_th_unmerge=1.4;                                                    
            distance_th=1.2;                                                            
            velocity_th=0.5;                                                            
            human_size_th=100;
            bus_size_th=2000;

            locksize_time=4;                                                            
            max_prediction_time=5;                                                      

            alpha_vel=0.8;                                                              
            alpha_size=0.5;                                                             
            alpha_dist=0;                                                               
            alpha_rand=0;                                                               
            max_prediction_time=0;                                                      

            alpha_dist=0;                                                               
            alpha_rand=0;                                                               

            alpha_vel=0;                                                              
            alpha_size=0;                                                             

            for i=1:num_tracker
                lastmatch{i}.pos=zeros(1,4);
                feature{i}.memwrite=0;
                feature{i}.testfeature0=zeros(1,128);
                feature{i}.testfeature1=zeros(1,128);
            end

            
            blank.locked=zeros(1,num_tracker);                                          
            blank.xvelocity=zeros(1,num_tracker);                                       
            blank.yvelocity=zeros(1,num_tracker);                                       
            blank.x=zeros(1,num_tracker);                                               
            blank.y=zeros(1,num_tracker);                                               
            blank.w=zeros(1,num_tracker);                                               
            blank.h=zeros(1,num_tracker);                                               
            blank.id=zeros(1,num_tracker);                                              
            blank.history=zeros(1,num_tracker);                                         
            blank.match=zeros(1,num_tracker);                                           
            blank.lockprediction=zeros(1,num_tracker);                                    
            blank.change=zeros(1,num_tracker);
            blank.dontmerge=zeros(1,num_tracker);
            blank.merge=zeros(1,num_tracker);
            blank.locksize=zeros(1,num_tracker);
            blank.refwidth=zeros(1,num_tracker);
            blank.refheight=zeros(1,num_tracker);
            blank.occlusion=zeros(1,num_tracker);
            blank.classify=zeros(1,num_tracker);
            blank.car=zeros(1,num_tracker);
            blank.bus=zeros(1,num_tracker);
            blank.bike=zeros(1,num_tracker);
            blank.truck=zeros(1,num_tracker);
            blank.human=zeros(1,num_tracker);
            blank.unknown=zeros(1,num_tracker);
            blank.occlusion_w=zeros(1,num_tracker);
            blank.occlusion_h=zeros(1,num_tracker);
            blank.entry_exit=zeros(1,num_tracker);

            
            D=zeros(num_tracker,num_tracker);
            V_D=zeros(num_tracker,num_tracker);
            Sim=zeros(num_tracker,num_tracker);

            past=blank;                                                                
            current=blank;
            track_id=0;            
                      
            tindex_count = 1;
            for jj = 1:length(TD.ts)
                t = TD.ts(jj);
                x = TD.x(jj);
                y = TD.y(jj);
                frame_temp(y,x) = 1;
                
                if (t-tstart)>obj.consts.window_period
                    
                    frame_ctr = frame_ctr+1;
                    i = frame_ctr;
                    bin_image{frame_ctr}.x = zeros(1,2*defnum_region);
                    bin_image{frame_ctr}.y = zeros(1,2*defnum_region);
                    bin_image{frame_ctr}.num_obj = 0;                  
                    
                    frame_temp2 = frame_temp;
                    frame_66ms = imresize(frame_temp2,[max_y,max_x]);
                   filt_frame=my_medfilt2(frame_66ms,...
                                  obj.consts.stride); 
                   bin_image{frame_ctr}.frame = filt_frame;
                             
                    if strcmp(RP_type, "histogram_based_RP")
                        [tracker_boxes, N_boxes] = ...
                                histogram_based_RPDVSslens(obj.consts,...
                                               filt_frame,t);
                    end
                    show_boxes5(tindex_count,filt_frame,0,tracker_boxes);
                    
                    if size(tracker_boxes,2) > defnum_region
                       end_num = defnum_region;
                       for start_num = 1:1:end_num
                           bin_image{frame_ctr}.x(2*(start_num-1)+1) = tracker_boxes(start_num).x;
                           bin_image{frame_ctr}.x(2*(start_num-1)+2) = tracker_boxes(start_num).x+tracker_boxes(start_num).xSize;
                           bin_image{frame_ctr}.y(2*(start_num-1)+1) = tracker_boxes(start_num).y;
                           bin_image{frame_ctr}.y(2*(start_num-1)+2) = tracker_boxes(start_num).y+tracker_boxes(start_num).ySize; 
                       end                   
                       bin_image{frame_ctr}.num_obj = end_num;
                    end
                    if size(tracker_boxes,2) <= defnum_region && size(tracker_boxes,2) ~= 0
                       end_num = size(tracker_boxes,2);
                       for start_num = 1:1:end_num
                           bin_image{frame_ctr}.x(2*(start_num-1)+1) = tracker_boxes(start_num).x;
                           bin_image{frame_ctr}.x(2*(start_num-1)+2) = tracker_boxes(start_num).x+tracker_boxes(start_num).xSize;
                           bin_image{frame_ctr}.y(2*(start_num-1)+1) = tracker_boxes(start_num).y;
                           bin_image{frame_ctr}.y(2*(start_num-1)+2) = tracker_boxes(start_num).y+tracker_boxes(start_num).ySize; 
                       end
                       bin_image{frame_ctr}.num_obj = end_num;
                    end              
                    
                    if (N_boxes > 0)                       
                           
                           frame_66ms = single(zeros(max_y, max_x));
                           tstart=TD.ts(jj);  
                           frame_temp = single(zeros(260,346));

                    end
                    if (N_boxes == 0) 
                        frame_66ms = single(zeros(max_y, max_x));
                        tstart=TD.ts(jj);  
                        frame_temp = single(zeros(260,346));
                    end
 
                    past=current;                                                           
                    current=blank;                                                          
                    empty_assigned=0;                                                       
                    ind_locked=find(past.locked==2);                                        
                    ind_tracking=find(past.locked==1);                                      
                    ind_empty=find(past.locked==0);                                         
                    matched=zeros(num_region,num_tracker);                                  
                    
                    num_obj = bin_image{i}.num_obj;                                         
                    if num_obj > num_tracker                                                          
                        num_obj = num_tracker;                                                        
                    end
                    for j=1:num_obj                                                         

                        R_P=[bin_image{i}.x(2*j-1) bin_image{i}.y(2*j-1) (bin_image{i}.x(2*j)-bin_image{i}.x(2*j-1)) (bin_image{i}.y(2*j)-bin_image{i}.y(2*j-1))]; 
                        indictor_RP = notemptyregion(R_P,bin_image{i}.frame);
                        if indictor_RP == 0
                           continue;
                        end          
                        if isempty(ind_locked)==0                                                                        

                            for k=1:length(ind_locked)                                                                   
                                q=ind_locked(k);                                                                         
                                
                                T_P=[past.x(q)+past.xvelocity(q) past.y(q)+past.yvelocity(q) past.w(q) past.h(q)];       
                                
                                overlap=rectint(T_P,R_P);                                                                
                                area_tracker=T_P(3)*T_P(4);                                                              
                                area_region=R_P(3)*R_P(4);                                                               
                                if overlap/area_tracker>overlap_th || overlap/area_region>overlap_th                     
                                    matched(j,q)=1;                                                                      
                                    if current.match(q)==0                                                               

                                        current.locked(q)=2;                                                             
                                        current.history(q)=past.history(q)+1;                                                             
                                        current.x(q)=(1-alpha_size)*R_P(1)+alpha_size*T_P(1);                            
                                        current.y(q)=(1-alpha_size)*R_P(2)+alpha_size*T_P(2);                            
                                        current.w(q)=(1-alpha_size)*R_P(3)+alpha_size*T_P(3);                            
                                        current.h(q)=(1-alpha_size)*R_P(4)+alpha_size*T_P(4);                            
                                        current.id(q)=past.id(q);                                                        
                                        current.match(q)=1;                                                              
                                        lastmatch{q}.pos=R_P;                                                            
                                                                              
                                        v1=R_P(1)- T_P(1);                                                              
                                        v2=R_P(1)+R_P(3) - T_P(1) - T_P(3);
                                        if abs(v1-past.xvelocity(q))<abs(v2-past.xvelocity(q))
                                            current.xvelocity(q)=v1;
                                        else
                                            current.xvelocity(q)=v2;
                                        end
                                        v1=R_P(2)- T_P(2);                                                                
                                        v2=R_P(2)+R_P(4) - T_P(2) - T_P(4);
                                        if abs(v1-past.yvelocity(q))<abs(v2-past.yvelocity(q))
                                            current.yvelocity(q)=v1;
                                        else
                                            current.yvelocity(q)=v2;
                                        end
                                        if past.xvelocity(q)~=0
                                            current.xvelocity(q)=(1-alpha_vel)*current.xvelocity(q)+alpha_vel*past.xvelocity(q);     
                                        end
                                        if past.yvelocity(q)~=0
                                            current.yvelocity(q)=(1-alpha_vel)*current.yvelocity(q)+alpha_vel*past.yvelocity(q);     
                                        end

                                    else                                                                                             
                                        
                                        distance=rectdist(R_P,lastmatch{q}.pos);                                                     
                                        
                                        current.dontmerge(q)=0;
                                        merge=rectmerge(R_P,[current.x(q) current.y(q) current.w(q) current.h(q)]);                  
                                        current.x(q)=(1-alpha_size)*merge(1)+alpha_size*T_P(1);                                      
                                        current.y(q)=(1-alpha_size)*merge(2)+alpha_size*T_P(2);                                      
                                        current.w(q)=(1-alpha_size)*merge(3)+alpha_size*T_P(3);                                      
                                        current.h(q)=(1-alpha_size)*merge(4)+alpha_size*T_P(4);                                      

                                        v1=merge(1)- T_P(1);                                                                          
                                        v2=merge(1)+merge(3) - T_P(1) - T_P(3);
                                        if abs(v1-past.xvelocity(q))<abs(v2-past.xvelocity(q))
                                            current.xvelocity(q)=v1;
                                        else
                                            current.xvelocity(q)=v2;
                                        end
                                        v1=merge(2)- T_P(2);                                                                          
                                        v2=merge(2)+merge(4) - T_P(2) - T_P(4);
                                        if abs(v1-past.yvelocity(q))<abs(v2-past.yvelocity(q))                                        
                                            current.yvelocity(q)=v1;
                                        else
                                            current.yvelocity(q)=v2;
                                        end
                                        if past.xvelocity(q)~=0
                                            current.xvelocity(q)=(1-alpha_vel)*current.xvelocity(q)+alpha_vel*past.xvelocity(q);     
                                        end
                                        if past.yvelocity(q)~=0
                                            current.yvelocity(q)=(1-alpha_vel)*current.yvelocity(q)+alpha_vel*past.yvelocity(q);     
                                        end
                                        
                                    end

                                end
                                indicator_current = notemptyregion2(current,bin_image{i}.frame,q);
                                if indicator_current == 0
                                   current.locked(q)=0;
                                   current.x(q)=0;
                                   current.y(q)=0;
                                   current.w(q)=0;
                                   current.h(q)=0;
                                   current.match(q)=0;
                                   current.xvelocity(q)=0;
                                   current.yvelocity(q)=0;
                                end
                            end


                        end

                        if sum(matched(j,:))>0
                            continue; 
                        end

                        if isempty(ind_tracking)==0                                                                                   
                            
                            for k=1:length(ind_tracking)
                                q=ind_tracking(k);                                                                                    
                                T_P=[past.x(q) past.y(q) past.w(q) past.h(q)]; 
                                   
                                overlap=rectint(T_P,R_P);
                                area_tracker=T_P(3)*T_P(4);
                                area_region=R_P(3)*R_P(4);
                                if overlap/area_tracker>overlap_th || overlap/area_region>overlap_th
                                    matched(j,q)=1; 
                                    if current.match(q)==0 

                                        current.locked(q)=2;                                                                          
                                        current.x(q)=(1-alpha_size)*R_P(1)+alpha_size*T_P(1);
                                        current.y(q)=(1-alpha_size)*R_P(2)+alpha_size*T_P(2);
                                        current.w(q)=(1-alpha_size)*R_P(3)+alpha_size*T_P(3);
                                        current.h(q)=(1-alpha_size)*R_P(4)+alpha_size*T_P(4);
                                        track_id=track_id+1;
                                        current.id(q)=track_id;
                                        current.match(q)=1;
                                       
                                        v1=R_P(1)- T_P(1);
                                        v2=R_P(1)+R_P(3) - T_P(1) - T_P(3);
                                        if abs(v1-past.xvelocity(q))<abs(v2-past.xvelocity(q))
                                            current.xvelocity(q)=v1;
                                        else
                                            current.xvelocity(q)=v2;
                                        end
                                        v1=R_P(2)- T_P(2);
                                        v2=R_P(2)+R_P(4) - T_P(2) - T_P(4);
                                        if abs(v1-past.yvelocity(q))<abs(v2-past.yvelocity(q))
                                            current.yvelocity(q)=v1;
                                        else
                                            current.yvelocity(q)=v2;
                                        end
                                        if past.xvelocity(q)~=0
                                            current.xvelocity(q)=(1-alpha_vel)*current.xvelocity(q)+alpha_vel*past.xvelocity(q);
                                        end
                                        if past.yvelocity(q)~=0
                                            current.yvelocity(q)=(1-alpha_vel)*current.yvelocity(q)+alpha_vel*past.yvelocity(q);
                                        end

                                    else                    
                                        merge=rectmerge(R_P,[current.x(q) current.y(q) current.w(q) current.h(q)]);
                                        current.x(q)=(1-alpha_size)*merge(1)+alpha_size*T_P(1);
                                        current.y(q)=(1-alpha_size)*merge(2)+alpha_size*T_P(2);
                                        current.w(q)=(1-alpha_size)*merge(3)+alpha_size*T_P(3);
                                        current.h(q)=(1-alpha_size)*merge(4)+alpha_size*T_P(4);

                                        v1=merge(1)- T_P(1);
                                        v2=merge(1)+merge(3) - T_P(1) - T_P(3);
                                        if abs(v1-past.xvelocity(q))<abs(v2-past.xvelocity(q))
                                            current.xvelocity(q)=v1;
                                        else
                                            current.xvelocity(q)=v2;
                                        end
                                        v1=merge(2)- T_P(2);
                                        v2=merge(2)+merge(4) - T_P(2) - T_P(4);
                                        if abs(v1-past.yvelocity(q))<abs(v2-past.yvelocity(q))
                                            current.yvelocity(q)=v1;
                                        else
                                            current.yvelocity(q)=v2;
                                        end
                                        if past.xvelocity(q)~=0
                                            current.xvelocity(q)=(1-alpha_vel)*current.xvelocity(q)+alpha_vel*past.xvelocity(q);
                                        end
                                        if past.yvelocity(q)~=0
                                            current.yvelocity(q)=(1-alpha_vel)*current.yvelocity(q)+alpha_vel*past.yvelocity(q);
                                        end
                                    end

                                end
                                indicator_current = notemptyregion2(current,bin_image{i}.frame,q);
                                if indicator_current == 0
                                   current.locked(q)=0;
                                   current.x(q)=0;
                                   current.y(q)=0;
                                   current.w(q)=0;
                                   current.h(q)=0;
                                   current.match(q)=0;
                                   current.xvelocity(q)=0;
                                   current.yvelocity(q)=0;
                                end
                            end

                        end
                        if sum(matched(j,:))>0
                            continue; 
                        end
                        if isempty(ind_empty)==0 && length(ind_empty)>empty_assigned                                                     
                         
                            empty_assigned=empty_assigned+1;
                            q=ind_empty(empty_assigned);                                                                                                                  
                            current.locked(q)=1;
                            current.x(q)=R_P(1);
                            current.y(q)=R_P(2);
                            current.w(q)=R_P(3);
                            current.h(q)=R_P(4);
                            current.match(q)=0;
                            current.xvelocity(q)=0;
                            current.yvelocity(q)=0;

                        else
                            disp('No available free trackers for new region proposal');
                        end

                    end



                    %% Cleaning up all the trackers %%%
                    temp.match=past.match-current.match;
                    if any(temp.match)
                        for k=1:length(temp.match)
                            if temp.match(k)==1
                                if past.lockprediction(k)<max_prediction_time+1 %continue predicting position assuming occlusion or insufficient information
                                    current.locked(k)=2;
                                    current.lockprediction(k)=past.lockprediction(k)+1;
                                    current.x(k)=past.x(k)+past.xvelocity(k);
                                    current.y(k)=past.y(k)+past.yvelocity(k);

                                    if current.x(k)>a_maxXSize || current.y(k)>a_maxYSize %object out of scene -- stop
                                        current.locked(k)=0;
                                        continue;
                                    end
                                    current.w(k)=past.w(k);
                                    current.h(k)=past.h(k);
                                    current.id(k)=past.id(k);
                                    current.match(k)=1;
                                    current.xvelocity(k)=past.xvelocity(k);
                                    current.yvelocity(k)=past.yvelocity(k);
                                    current.history(k)=past.history(k)+1;
                                else 
                                    
                                    current.locked(k)=0;
                                end
                            end
                        end
                    end

                    track{i}=current;
                    tindex_count = tindex_count+1;
                end
                
            end
            obj.max_track_ID_found = max_trk_i;    
        end   
        
        
 
        
        
        
        
        
        
        
        function obj = show_box_with_accuracy_demo(obj,TD,...
                                annots,RP_type,xem)
            unique_track_id = unique(annots.track_num);
            obj.files_explored = obj.files_explored + 1; 
            N_class = 4;
            for track_idx = 1:length(unique_track_id)
%             for track_idx = 1:3
                if (mod(track_idx,floor(length(unique_track_id)/10))==0)
                    disp(['done another 10% ' num2str(track_idx) '/'...
                    num2str(length(unique_track_id))]);
                end

                track_num = unique_track_id(track_idx);
                track_idx_prv=max(track_idx-1,1);
                track_idx_nxt=min(track_idx+1,length(unique_track_id));
                track_num_prv=unique_track_id(track_idx_prv);
                track_num_nxt=unique_track_id(track_idx_nxt);
                track = getTrackVec_fixed_interval_wo_occl(TD,...
                        annots,track_num,...
                        track_num_prv,track_num_nxt,...
                        obj.consts.window_period,...
                        obj.consts.wo_entry_exit);            
                track = RemoveNulls(track, track.isPartial);
                classID  = track.meta.class;
                if (~isempty(track.ts) &&  (classID>0) && (classID<7) ...
                    && (classID~=4))
                    if (classID==3)
                        track.meta.class = 1;
                    end 
                    if (classID>4)
                        track.meta.class = classID - 2;
                    end
                    obj.track_for_demo(track,...
                                RP_type,xem,...
                                N_class);
                end
                
            end
        end

        
        function obj = show_box_with_accuracy_demo2(obj,TD,...
                                annots,RP_type,~)
            unique_track_id = unique(annots.track_num);
            obj.files_explored = obj.files_explored + 1; 
            N_class = 4;
            for track_idx = 1:length(unique_track_id)
%             for track_idx = 1:3
                if (mod(track_idx,floor(length(unique_track_id)/10))==0)
                    disp(['done another 10% ' num2str(track_idx) '/'...
                    num2str(length(unique_track_id))]);
                end

                track_num = unique_track_id(track_idx);
                track_idx_prv=max(track_idx-1,1);
                track_idx_nxt=min(track_idx+1,length(unique_track_id));
                track_num_prv=unique_track_id(track_idx_prv);
                track_num_nxt=unique_track_id(track_idx_nxt);
                track = getTrackVec_fixed_interval_wo_occl(TD,...
                        annots,track_num,...
                        track_num_prv,track_num_nxt,...
                        obj.consts.window_period,...
                        obj.consts.wo_entry_exit);            
                track = RemoveNulls(track, track.isPartial);
                classID  = track.meta.class;
%                 if (~isempty(track.ts) &&  (classID>0) && (classID<7) ...
%                     && (classID~=4))
%                     if (classID==3)
%                         track.meta.class = 1;
%                     end 
%                     if (classID>4)
%                         track.meta.class = classID - 2;
%                     end
%                     obj.track_for_demo(track,...
%                                 RP_type,xem,...
%                                 N_class);
%                 end
                
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        function obj = explore_video_usign_GT(obj,TD,annots)
            unique_track_id = unique(annots.track_num);
            obj.files_explored = obj.files_explored + 1; 
            for track_idx = 1:length(unique_track_id)
%             for track_idx = 1:3
                if (mod(track_idx,floor(length(unique_track_id)/10))==0)
                    disp(['done another 10% ' num2str(track_idx) '/'...
                    num2str(length(unique_track_id))]);
                end

                track_num = unique_track_id(track_idx);
                track_idx_prv=max(track_idx-1,1);
                track_idx_nxt=min(track_idx+1,length(unique_track_id));
                track_num_prv=unique_track_id(track_idx_prv);
                track_num_nxt=unique_track_id(track_idx_nxt);
                track = getTrackVec_fixed_interval_wo_occl(TD,...
                        annots,track_num,...
                        track_num_prv,track_num_nxt,...
                        obj.consts.window_period,...
                        obj.consts.wo_entry_exit);            
                track = RemoveNulls(track, track.isPartial);
                [track_hist_feat, track_hist_feat_sizes,...
                track_aug_hist_feat,...
                track_length,no_aug_ims,...
                istrackempty]  = obj.explore_track(track);
                if (~istrackempty) && (track.meta.class>0) && ...
                    (track.meta.class<=obj.consts.no_classes)
                    obj.current_track_ID = obj.current_track_ID + 1;
                    % all_img_filt = [all_img_filt, uint16([ones(1,...
                    %      size(img_filt,2)) .* cur_class ; img_filt])];
                    % all_img = [all_img, uint16([ones(1, size(img,2)) .* ...
                      % cur_class ; img])];
                    obj.add_features(track_hist_feat,...
                        track_hist_feat_sizes,...
                        ones(1,track_length)*track.meta.class,...
                        track_length,zeros(1,track_length),obj.files_explored);
                    obj.all_examples = obj.all_examples + track_length;
                    obj.examples_by_class(track.meta.class) = ...
                        obj.examples_by_class(track.meta.class) + track_length;
                    obj.tracks_by_class(track.meta.class) = ...
                        obj.tracks_by_class(track.meta.class) + 1;
                    if no_aug_ims>0
                       aug_hst_sizes  = zeros(2,no_aug_ims);
                       obj.add_features(track_aug_hist_feat,...
                            aug_hst_sizes,...
                            ones(1,no_aug_ims)*track.meta.class,...
                            no_aug_ims,ones(1,no_aug_ims),obj.files_explored);
                       obj.all_examples = obj.all_examples + no_aug_ims;
                       obj.examples_by_class(track.meta.class) = ...
                           obj.examples_by_class(track.meta.class) + no_aug_ims;
                    end
                end
            % all_hist_features_vel_dir = [all_hist_features_vel_dir,...
            %                              ones(1,size(hist_features,2))*vel_dir];
            end
        end
    end
    
    methods(Access = private)
        function add_features(obj,...
                 track_hist_feat,...
                 track_hist_feat_sizes,...
                 class_ID,track_length,...
                 augmented,fileID)
            % Not sure whether concatenation will be correct
            obj.hist_features = [obj.hist_features, ...
                                track_hist_feat];
            obj.meta_data.add_meta_data(...
                ones(1,track_length)*obj.current_track_ID,...
                track_hist_feat_sizes(1,:),...
                track_hist_feat_sizes(2,:),class_ID,...
                augmented,track_length,fileID);

        end
        function track_for_demo(obj,TD,...
                                RP_type,xem,...
                                N_class)
            max_y = obj.consts.maxY;
            max_x = obj.consts.maxX;
            tstart=TD.ts(1);
            max_hist_size = obj.consts.max_hist_size;
            frame_66ms = single(zeros(max_y, max_x));
            classID  = TD.meta.class;
            colors = {'red','green','blue','cyan',...
                        'magenta','yellow','red','green'};
            cats = {'CarVan', 'Bus','Bike','Truck'};
            X = categorical(cats);
            X = reordercats(X,cats);
            % writerObj = VideoWriter('file_video.avi', 'Uncompressed AVI');
            %  writerObj.FrameRate = 24;
            % open(writerObj);
            frames_encounterd = 0;
            voted_arr = zeros(1,4);
            for jj = 1:length(TD.ts)
                t = TD.ts(jj);
                x = TD.x(jj);
                y = TD.y(jj);
                frame_66ms(y,x) = 1;
                
                if (t-tstart)>obj.consts.window_period
                    
                    filt_frame=...
                                my_medfilt2(frame_66ms,...
                                obj.consts.stride); 
                    
                    if strcmp(RP_type, "histogram_based_RP")
                        [tracker_boxes, N_boxes] = ...
                                histogram_based_RP(obj.consts,...
                                               filt_frame,t);
                    elseif strcmp(RP_type, "CCA_based_RP")
                        [tracker_boxes, N_boxes] = ...
                                Single_Pass_CCL(obj.consts,...
                                               filt_frame,t);
                    end

                    frame_66ms_copy = frame_66ms;
                    frame_66ms_copy(frame_66ms_copy==0) = 255;
                    frame_66ms_copy(frame_66ms_copy==1) = 0;
                    image2 = zeros(max_y,...
                                max_x,3);
        
                    for kk = 1:3
                        image2(:,:,kk) = frame_66ms_copy;
                    end
        
                    image2 = uint8(image2);
                    
                    max_No = min(8,N_boxes);
                    votes = [];
                    figure(1);
                    for kk = 1:max_No
                        x = tracker_boxes(kk).x;
                        y = tracker_boxes(kk).y;
                        xSize = tracker_boxes(kk).xSize;
                        ySize = tracker_boxes(kk).ySize;
                        A= [x,y,xSize,ySize];
                        image2 = insertShape(image2,'rectangle',...
                            A,'Color',colors(kk));
                        this_box = tracker_boxes(kk);
                        object = frame(this_box.y:(this_box.y + this_box.ySize-1),...
                               this_box.x:(this_box.x + this_box.xSize-1)); 
                        [hist_features, ...
                             hist_features_sizes]=...
                                  hist_FE(obj.consts.hist_consts,...
                                   max_hist_size,object);
                        hist_features  = hist_features - 128;
                        probabs = pass_hist_from_chip(xem,...
                            hist_features,...
                            classID,N_class);
                        subplot(3,4,4+kk);
                        bar(X,probabs,colors(kk));
                        title('Region Proposal %d',kk);
                        [Maxprob,class_mine] = max(probabs);
                        votes = [votes,class_mine];
                    end
                    voted_class = mode(votes);
                    voted_arr(voted_class) = voted_arr(voted_class) + 1;
                    subplot(3,4,[2,3]);
                    imshow(image2);
                    title('Current Image');
                    axis on;

                    subplot(3,4,[10,11]);
                    bar(X,voted_arr,colors(kk));
                    title('Voting Accuracy');
                    pause;
                    frame_66ms = single(zeros(max_y, max_x));
                    tstart=TD.ts(jj);

                end
                
            end
            close(writerObj);
        end
        function add_features_diff_trks(obj,...
                 track_hist_feat,...
                 track_hist_feat_sizes,...
                 class_ID,trk_IDs,track_length,...
                 augmented,fileID)
            % Not sure whether concatenation will be correct
            obj.hist_features = [obj.hist_features, ...
                                track_hist_feat];
            obj.meta_data.add_meta_data(...
                trk_IDs,...
                track_hist_feat_sizes(1,:),...
                track_hist_feat_sizes(2,:),class_ID,...
                augmented,track_length,fileID);

        end
        function [hist_features, hist_features_sizes,...
                aug_hist_features,aug_cpp_classes,...
                aug_cpp_trackIDs,aug_imgs_ctr] = ...
                                            explore_boxes(obj,frame,...
                                            boxes,cpp_classes,...
                                            cpp_trackIDs,N_boxes)
            max_hist_size = obj.consts.max_hist_size;

            hist_features = single(zeros(...
                                2*max_hist_size, N_boxes));
            hist_features_sizes = zeros(2, N_boxes);
            aug_hist_features = [];
            aug_cpp_classes = [];
            aug_cpp_trackIDs = [];
            aug_imgs_ctr = 0;
            for frame_ctr = 1:N_boxes
                this_box = boxes(frame_ctr);
                object = frame(this_box.y:(this_box.y + this_box.ySize-1),...
                               this_box.x:(this_box.x + this_box.xSize-1)); 
                [hist_features(:, frame_ctr), ...
                             hist_features_sizes(1:2,frame_ctr)]=...
                                  hist_FE(obj.consts.hist_consts,...
                                   max_hist_size,object);
                min_size = min(hist_features_sizes(:,frame_ctr));
                if (ismember(cpp_classes(frame_ctr),...
                    obj.consts.classes_to_augment) && (min_size>2))
                    [aug_hist_feat, n_aug_img] = augment_image(...
                          obj.consts.hist_consts,...
                          object,max_hist_size,...
                          "hist_FE");
                    aug_hist_features = [aug_hist_features,...
                        aug_hist_feat];
                    aug_cpp_classes = [aug_cpp_classes,...
                                        ones(1,n_aug_img)*...
                                        cpp_classes(frame_ctr)];
                    aug_cpp_trackIDs = [aug_cpp_trackIDs,...
                                        ones(1,n_aug_img)*...
                                        cpp_trackIDs(frame_ctr)];                    
                    aug_imgs_ctr = aug_imgs_ctr + n_aug_img;
                end
            end
        end

        function [hist_features, hist_features_sizes,...
                  aug_hist_features,track_length,...
                  aug_imgs_ctr,istrackempty] = explore_track(obj,track)
            hist_features = [];
            hist_features_sizes = [];
            aug_hist_features = [];
            if ~isempty(track.ts)
                istrackempty = false;
                max_hist_size = obj.consts.max_hist_size;
                max_y = obj.consts.maxY;
                max_x = obj.consts.maxX;
                track.ts=track.ts-track.ts(1)+1;
% 
                num_hist_features = floor(max(track.ts) /...
                    obj.consts.window_period);

                hist_features = single(zeros(...
                                2*max_hist_size, num_hist_features));
                % What is the purpose of this img variable?
                % img = single(zeros(max_hist_size*max_hist_size, num_hist_features));

                hist_features_sizes = zeros(2, num_hist_features);
% 
                frame_ctr=1;
                aug_imgs_ctr = 0;
                tstart=track.ts(1);
                prev_jj=1;
                frame_33ms = single(zeros(max_y, max_x));
                frame_33ms_prev=frame_33ms;
                for jj=1:length(track.ts)

                    frame_33ms(track.y(jj), track.x(jj)) = 1;
                    t=track.ts(jj);
                    if t-tstart>obj.consts.window_period 

                        consider_Y=find(sum(frame_33ms,1)>0);
                        consider_X=find(sum(frame_33ms,2)>0);
% 
%                     %% Applies the median filter over the frame so as to 
%                     % remove the static noise. And also the size of the 
%                     % image can be changed.
                        if (length(consider_X)>=2) && (length(consider_Y)>=2)
                            frame_33ms(consider_X(1):consider_X(end),...
                                consider_Y(1):consider_Y(end))=...
                                    my_medfilt2(frame_33ms(consider_X(1):...
                                    consider_X(end),...
                                    consider_Y(1):consider_Y(end)),...
                                    obj.consts.stride);
                        end
                        consider_Y=find(sum(frame_33ms,1)>0);
                        consider_X=find(sum(frame_33ms,2)>0);
 
                        if (length(consider_X)>=2) && (length(consider_Y)>=2)
                            frame=frame_33ms(consider_X(1):consider_X(end),...
                                consider_Y(1):consider_Y(end));
%                             figure(1),imshow(frame);
%                             hold on;
%                             rectangle('Position',[consider_Y(1),consider_X(1),...
%                                 consider_Y(end)-consider_Y(1)+1,...
%                                 consider_X(end)-consider_X(1)+1],'Edgecolor',...
%                                 'g');
                            [hist_features(:, frame_ctr) ...
                             hist_features_sizes(1:2,frame_ctr)]=...
                                  hist_FE(obj.consts.hist_consts,...
                                   max_hist_size,frame);
                            if ismember(track.meta.class,...
                                obj.consts.classes_to_augment) 
                                [aug_hist_feat, n_aug_img] = augment_image(...
                                      obj.consts.hist_consts,...
                                      frame,max_hist_size,...
                                      "hist_FE");
                                aug_hist_features = [aug_hist_features,...
                                    aug_hist_feat];
                                aug_imgs_ctr = aug_imgs_ctr + n_aug_img;
                            end
                            frame_ctr=frame_ctr+1;
                        end
% 
                    frame_33ms_prev=frame_33ms;
                    frame_33ms = single(zeros(max_y, max_x));
                    tstart=track.ts(jj);
                    prev_jj=jj;
                    end
                    
                end
                track_length = frame_ctr - 1;
                hist_features = hist_features(:,1:track_length);
                hist_features_sizes = hist_features_sizes(:,1:track_length);
            else
                istrackempty = true; 
                track_length = 0; 
                aug_imgs_ctr = 0;  
            end

        end
    end
    
    methods (Static)
        function obj = init_database
            obj.meta_data = meta_data_class.init_meta_data;
            obj.hist_features = [];
            obj.img_features  = [];
        end
    end
end
