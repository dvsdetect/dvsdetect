function hist_scaled = get_first_nhist(histogram, max_size)
	hist_scaled = zeros(1,max_size);
	if length(histogram)>=max_size
		hist_scaled(:) = histogram(1:max_size);
	else
		hist_scaled(1:length(histogram)) = histogram;	
	end
end
