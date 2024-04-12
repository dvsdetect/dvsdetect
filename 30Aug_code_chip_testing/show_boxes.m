function  RGB = show_boxes(filt_frame,gt_boxes,...
                                        tracker_boxes,...
                                        tracker_boxes2)

    filt_frame(filt_frame==0) = 255;
    filt_frame(filt_frame==1) = 0;
    [H,W] = size(filt_frame); 
    if (nargin ==3)
    figure(1), imshow(filt_frame);
    axis on;
	hold on;
% 	for ii = 1:length(gt_boxes)
% 		x = gt_boxes(ii).x;
% 		y = gt_boxes(ii).y;
% 		xSize = gt_boxes(ii).xSize;
% 		ySize = gt_boxes(ii).ySize;
% 		y = y - ySize + 1;
% 		rectangle('Position',[x,y,xSize,ySize],'Edgecolor', 'g');
% 
% 	end

	for ii = 1:length(tracker_boxes)
		x = tracker_boxes(ii).x;
		y = tracker_boxes(ii).y;
		xSize = tracker_boxes(ii).xSize;
		ySize = tracker_boxes(ii).ySize;
		% y = y - ySize + 1;
		rectangle('Position',[x,y,xSize,ySize],'Edgecolor', 'r');
		
    end
    elseif (nargin==4)
        image2 = zeros(2*H,...
                        W,3);
        for ii= 1:2
            for jj = 1:3
                
            image2((ii-1)*H +1:ii*H,:,jj) = filt_frame;
            end
        end
        image2(H-1:H+1,:,1) = 0;
        image2(H-1:H+1,:,2) = 0;
        image2(H-1:H+1,:,3) = 255;
        image2 = uint8(image2);
        A =[];
        
        for ii = 1:length(tracker_boxes)
            x = tracker_boxes(ii).x;
            y = tracker_boxes(ii).y;
            xSize = tracker_boxes(ii).xSize;
            ySize = tracker_boxes(ii).ySize;
            % y = y - ySize + 1;
            A= [A;[x,y,xSize,ySize]];
            
        end
        for ii = 1:length(tracker_boxes2)
            x = tracker_boxes2(ii).x;
            y = tracker_boxes2(ii).y;
            xSize = tracker_boxes2(ii).xSize;
            ySize = tracker_boxes2(ii).ySize;
            % y = y - ySize + 1;
            A= [A;[x,y+H,xSize,ySize]];
            
        end
        RGB = insertShape(image2,'rectangle',...
                A,'Color','red');
        
        imshow(RGB)
    end
end
