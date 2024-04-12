function [score1, score2,score3]= compute_overlap_score2(gt_box, eval_box)

gt_box.x=double(gt_box.x);
eval_box.x=double(eval_box.x);
gt_box.xSize=double(gt_box.xSize);
eval_box.xSize=double(eval_box.xSize);
gt_box.y=double(gt_box.y);
eval_box.y=double(eval_box.y);
gt_box.ySize=double(gt_box.ySize);
eval_box.ySize=double(eval_box.ySize);
% compute the vertices of the overlapping area of the rectangles
left    = max(gt_box.x , eval_box.x );
right   = min(gt_box.x + gt_box.xSize, eval_box.x + eval_box.xSize);
top     = max(gt_box.y , eval_box.y);
bottom  = min(gt_box.y + gt_box.ySize, eval_box.y + eval_box.ySize);

% compute overlapping area (intersection)
if (left<right && bottom > top)
    intersection = (right-left)*(bottom-top);
else
    intersection = 0;
end

%compute area of each rectangle
gt_area = gt_box.xSize*gt_box.ySize;
eval_area = eval_box.xSize*eval_box.ySize;

% compute the intersection over the union
union = gt_area+eval_area-intersection;

score1 = intersection/union;

if (union== gt_area || union == eval_area)
    
    score2 = 1;
else
    
    score2 = intersection/gt_area;
    %score2 = intersection/eval_area;
end

if (union== gt_area || union == eval_area)
    
    score3 = 1;
else
    
    score3 = intersection/eval_area;
    %score2 = intersection/eval_area;
end

% if (score2<score1)
%     disp([score1 score2]);
%     disp([gt_area intersection union]);
% end


end

