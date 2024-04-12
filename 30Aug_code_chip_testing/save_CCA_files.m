function save_CCA_files(data1,data2)
    data1(data1==0) = 255;
    data1(data1==1) = 0;
    
    imshow(data1);
    axis on
    
    data2(data2==0) = 255;
    data2(data2==1) = 0;
    imshow(data2);
    axis on
    
    data3 = ones(size(data2,1),size(data2,2),3);
    data3(:,:,1) = ones(size(data2))*255;
    data3(:,:,2) = data2;
    data3(:,:,3) = data2;
    imshow(data3);
    axis on
end