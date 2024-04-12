function annotation = read_annotation(filename)
fileID = fopen(filename, 'r');

fscanf(fileID,'#This text file contains annotation data for recordings in: ');
annotation.meta.source = fgetl(fileID); %for the newline character etc
annotation.meta.source = strrep(annotation.meta.source, '\', filesep);
annotation.meta.source = strrep(annotation.meta.source, '/', filesep);

fscanf(fileID, '#The corresponding picture of the recording site is at: ');
annotation.meta.picture = fgetl(fileID); % for the newline character etc
annotation.meta.picture = strrep(annotation.meta.picture, '\', filesep);
annotation.meta.picture = strrep(annotation.meta.picture, '/', filesep);

fscanf(fileID, '#The annotation file is stored at: ');
annotation.meta.annotation_file = fgetl(fileID);
annotation.meta.annotation_file = strrep(annotation.meta.annotation_file, '\', filesep);
annotation.meta.annotation_file = strrep(annotation.meta.annotation_file, '/', filesep);

fscanf(fileID,'#Comments: ');
annotation.meta.comments = fgetl(fileID); %for the newline character etc

fscanf(fileID,'#The recordings are annotated by: ');
annotation.meta.annotator = fgetl(fileID);

annotation.meta.sensor_dimension = flipud(fscanf(fileID,'#Sensor Dimensions- Height = %d Pixels Width = %d Pixels'));
fgets(fileID); %for the newline character etc
fscanf(fileID, '#LEGEND: ');
annotation.meta.legend = fgetl(fileID);
fgets(fileID); %for the newline character etc

A = fscanf(fileID,'%u, %d, %d, %u, %u, %u, %u,\n', [7,inf])'; %7 columns, all rows
annotation.ts = A(:,1);
annotation.x = A(:,2);
annotation.y = A(:,3);
annotation.xSize = A(:,4);
annotation.ySize = A(:,5);
annotation.track_num = A(:,6);
annotation.class = A(:,7);

fclose(fileID);
end

