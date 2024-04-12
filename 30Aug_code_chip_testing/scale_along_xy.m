function img_scaled = scale_along_xy(data,sx,sy)
    
    [height,width] = size(data);
    img_scaled = zeros(height/sy,width/sx);
    row_cn = 0;
    for row = 1:sy:height
        row_cn = row_cn + 1;
        col_cn = 0;
        for col = 1:sx:width
            col_cn = col_cn + 1;
            temp_ans = data(row:(row+sy-1),col:(col + sx -1));
            if (sum(temp_ans(:))>0)
                
                img_scaled(row_cn, col_cn) = 1;
            end
        end
    end
end

