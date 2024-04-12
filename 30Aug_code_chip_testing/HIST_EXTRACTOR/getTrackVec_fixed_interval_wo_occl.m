function track = getTrackVec_fixed_interval_wo_occl(TD, annotation,...
                 track_num,track_num_prv,track_num_nxt,fixed_interval,wo_entry_exit)
%function track = getTrackVec_fixed_interval(TD, annotation, track_num,fixed_interval, skip, skipby)
%function track = getTrack(TD, annotation, track_num)
% gets the tracker number "track_num" from the events passed
% Does not yet support dynamic sizing

% Removes the annotation except the one with the trank_num equal to track_num_*
annotation_prv = RemoveNulls(annotation, annotation.track_num ~=track_num_prv);
annotation_nxt = RemoveNulls(annotation, annotation.track_num ~=track_num_nxt);
annotation = RemoveNulls(annotation, annotation.track_num ~=track_num);

% sorts the inputEvents stream by temporal order using the 'ts' field or
% 'Timestamp' field.
annotation = SortOrder(annotation);

TD = RemoveNulls(TD, TD.ts<=annotation.ts(1)|TD.ts>=annotation.ts(end));
TD = SortOrder(TD);


occl_min=annotation.ts(1);
occl_max=annotation.ts(end);


%% Couldn't understand this part? 
if (wo_entry_exit==1)
    if (length(annotation.ts)>=4)
        occl_min=max(occl_min,annotation.ts(2));
        occl_max=min(occl_max,annotation.ts(end-1));
    end
end

max_x = max(TD.x);
max_y = max(TD.y);

TD.cenX = zeros(size(TD.ts)); % X-center of object with dynamic sizing
TD.cenY = zeros(size(TD.ts)); % Y-center of object with dynamic sizing
TD.xSize = zeros(size(TD.ts),'double'); % object size in pixels
TD.ySize = zeros(size(TD.ts),'double'); % object size in pixels
TD.isPartial = zeros(size(TD.ts),'uint8'); % 0 1 logical

TDinvalid = ones(size(TD.ts));
evt_start_idx = 1;

%% What's happening inside this for loop ??
for an_idx = 1:length(annotation.ts)-1
    evt_last_idx = find(TD.ts < annotation.ts(an_idx+1) ,1,'last');
    
    if (evt_last_idx>=evt_start_idx)
        %if ((skip==0) | (mod(an_idx,skipby)==0))
        
        TD_update_fixed_interval = annotation.ts(an_idx)+fixed_interval*...
                        floor((TD.ts(evt_start_idx:evt_last_idx)-annotation.ts(an_idx))/fixed_interval);
        time_factor = (TD_update_fixed_interval - annotation.ts(an_idx))/...
                     (annotation.ts(an_idx+1)- annotation.ts(an_idx)); %interpolate at fixed interval steps
        if (fixed_interval>10e6)
            time_factor(1:end)=0;   %don't interpolate at all;
        end
        posX = annotation.x(an_idx) + time_factor .*...
                 (annotation.x(an_idx + 1) - annotation.x(an_idx));
        posY = annotation.y(an_idx) + time_factor .*...
                 (annotation.y(an_idx + 1) - annotation.y(an_idx));
        xSize = annotation.xSize(an_idx) + time_factor .*...
                 (annotation.xSize(an_idx+1) - annotation.xSize(an_idx));
        ySize = annotation.ySize(an_idx) + time_factor .*...
                 (annotation.ySize(an_idx+1) - annotation.ySize(an_idx));
        
        inBox = and(abs(TD.x(evt_start_idx:evt_last_idx)-posX) < xSize/2, ...
            abs(TD.y(evt_start_idx:evt_last_idx)-posY) < ySize/2);
        TDinvalid(evt_start_idx:evt_last_idx) = not(inBox);
        
        left_boundary_x = fix(posX - xSize/2);
        left_boundary_y = fix(posY - ySize/2);
        right_boundary_x = fix(posX + xSize/2);
        right_boundary_y = fix(posY + ySize/2);
        
        isPartial = ...
            or(or(left_boundary_x < 1 ,...
            left_boundary_y < 1), ...
            or(right_boundary_x > max_x,...
            right_boundary_y > max_y));
        
        TD.isPartial(evt_start_idx:evt_last_idx) = isPartial;
        
        left_boundary_x = max(left_boundary_x,1);
        left_boundary_y = max(left_boundary_y,1);
        right_boundary_x = min(right_boundary_x,max_x);
        right_boundary_y = min(right_boundary_y,max_y);
        xSizeP = fix(right_boundary_x - left_boundary_x);
        ySizeP = fix(right_boundary_y - left_boundary_y);
        TD.xSize(evt_start_idx:evt_last_idx) = xSizeP .*...
                                         isPartial + fix(xSize) .* not(isPartial);
        TD.ySize(evt_start_idx:evt_last_idx) = ySizeP .*...
                                     isPartial + fix(ySize) .* not(isPartial);
        TD.cenX(evt_start_idx:evt_last_idx) = fix(left_boundary_x + xSizeP/2) .*...
                                         isPartial + fix(posX) .* not(isPartial);
        TD.cenY(evt_start_idx:evt_last_idx) = fix(left_boundary_y + ySizeP/2) .*...
                                             isPartial + fix(posY) .* not(isPartial);

    end
    
    evt_start_idx = evt_last_idx + 1;
end

TDinvalid(TD.ts<=occl_min)=1;
TDinvalid(TD.ts>=occl_max)=1;

track = RemoveNulls(TD, TDinvalid);
track = SortOrder(track);
track.meta.class = annotation.class(1);
track.meta.xSize = annotation.xSize(1);
track.meta.ySize = annotation.ySize(1);