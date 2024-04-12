function indictor_current = notemptyregion2(tracker_boxes,test_frame,q)
indictor_current = 0;
x = tracker_boxes.x(q);
if x<=0
   x=1;
end
y = tracker_boxes.y(q);  
if y<=0
   y=1;
end
xSize = tracker_boxes.w(q);
if xSize<=0
   xSize=1;
end
ySize = tracker_boxes.h(q);
if ySize<=0
    ySize=1;
end
if x+xSize>346
   xend = 346;
else
   xend = x+xSize;
end
if y+ySize>260
   yend = 260;
else
   yend = y+ySize;
end   
if sum(sum(test_frame(y:yend,x:xend)))>1
   indictor_current = 1;
end
end