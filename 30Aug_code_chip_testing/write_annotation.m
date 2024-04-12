classdef write_annotation < handle
    properties 
        consts
        fileOpen
        fileID 
        
    end
    methods
	    function obj = write_annotation(consts,filename)
	        if ~exist(consts.database_save_path,...
                    'dir')
                mkdir(consts.database_save_path);
            end
            formatOut = 'yyyymmdd';
            dt = datestr(now,formatOut);
            save_filename = [consts.database_save_path,'/',...
                            'annotations_',filename,'_',dt,'.txt'];
            fileID = fopen(save_filename, 'w');
			ini_comments = ['#This text file contains Kalman filter annotation data for recordings in: Recording\\',filename,'.bin\n', ...
			'#The corresponding picture of the recording site is at: Image\\ENG.png\n', ...
			'#The annotation file is stored at: Annotation\\',filename,'.txt\n', ...
			'#Comments: 12mm lens\n', ...
			'#The recordings are annotated by: ShiHao\n', ...
			'#Sensor Dimensions- Height = 180 Pixels Width = 240 Pixels\n', ...
			'#LEGEND: 1-Car, 2-Bus, 3-Van, 4-Pedestrian, 5-Bike, 6-Truck, 7-Unknown\n', ...
			'#Time(us),x-Location,y-Location,x-size,y-size,track-num,class\n'];
			fprintf(fileID, ini_comments);
			obj.consts = consts;
			obj.fileOpen = fileID;
			obj.fileID = filename;
		end

		function add_annotations(obj,tracker_boxes,...
							N_boxes,cpp_classes,cpp_trackIDs)
			for index = 1:N_boxes
				ts = tracker_boxes(index).ts;
				x = tracker_boxes(index).x;
				y = tracker_boxes(index).y;
				w = tracker_boxes(index).xSize;
				h = tracker_boxes(index).ySize;
				x = round(x + w/2);
				y = round(y + h/2);
				trackID = cpp_trackIDs(index);
				if (trackID<0)
					trackID = 0;
				end
				class = cpp_classes(index);
				fprintf(obj.fileOpen, '%ld, %d, %d, %d, %u, %u, %u,\n',...
                        ts, x, y, w, h, ...
					 trackID, class);
			end
		end

		function close(obj)
			fclose(obj.fileOpen);
		end
		
	end
end
