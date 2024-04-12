function  RGB = show_boxes5(t,filt_frame,gt_boxes,...
                                        tracker_boxes,...
                                        tracker_boxes2)

    filt_frame(filt_frame==0) = 255;
    filt_frame(filt_frame==1) = 0;
    [H,W] = size(filt_frame); 
    if (nargin ==4)
    fig = figure(1);
    figure(1), imshow(filt_frame);

	hold on;


	for ii = 1:length(tracker_boxes)
		x = tracker_boxes(ii).x;
		y = tracker_boxes(ii).y;
		xSize = tracker_boxes(ii).xSize;
		ySize = tracker_boxes(ii).ySize;
		rectangle('Position',[x,y,xSize,ySize],'Edgecolor', 'r','Linewidth',2);
    end
    set(gca,'position',[0 0 1 1]);
    axis normal
 
    elseif (nargin==5)
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
