function save_frame_with_hists(com_path,filt_bin_image,xhist,yhist)
	filt_img = filt_bin_image;
%     figure(1),imshow(filt_img);
	filt_img(filt_img==0) = 255;
    filt_img(filt_img==1) = 0;
    
	figure(1),imshow(filt_img);
	h = gca;
    h.Visible = 'On';
    imwrite(filt_img,[com_path,'/default_img.png']);

	figure(1)
    
    plot(1:length(xhist),xhist,'LineWidth',2)
    axis([1 length(xhist) 0 30])
    filename = sprintf([com_path,'/x_histogram_file.png']);
    saveas(gcf,filename)

	figure(1)
    
    plot(1:length(yhist),yhist,'LineWidth',2)
    axis([1 length(yhist) 0 30])
    filename = sprintf([com_path,'/y_histogram_file.png']);
    saveas(gcf,filename)    
end
