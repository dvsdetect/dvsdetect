function new_img_min=my_medfilt3(frame)
expand_size = 1;
img = wextend('2D','zpd',frame,expand_size);

[M,N] = size(frame);

muban_size = 3;
muban = ones(muban_size,muban_size);

for i = 1:1:M
    for j = 1:N
        mat = img(i:i+muban_size-1,j:j+muban_size-1);
        mat = mat(:);
        new_img_min(i,j) = min(mat);
    end
end