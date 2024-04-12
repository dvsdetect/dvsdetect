function a=NotToBeEliminated(annot_filename)

% load('/datasets/NEUPROWLER/Traffic_Scenes_NUS_NEW/TrueNorth_filtered/file_names_elimination.mat');
% 
% a=1;
% % file_names_train{1,1}
% % annot_filename
% for i=1:length(file_names_test)
%     if (strcmp(file_names_test{1,i},annot_filename))
%         a=0;
%     end
% end
% for i=1:length(file_names_train)
%     if (strcmp(file_names_train{1,i},annot_filename))
%         a=0;
%     end
% end
% for i=1:length(file_names_validate)
%     if (strcmp(file_names_validate{1,i},annot_filename))
%         a=0;
%     end
% end


filename = '/code/aakash/NEUPROWLER_rep_aakash/file_names_elimination.csv';
delimiter = '';
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
filenameselimination = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID ans;

a=1;
% filenameselimination{1,1}
% annot_filename
% filenameselimination{1,1}=annot_filename;
% filenameselimination{1,1}
for i=1:length(filenameselimination)
    if (strcmp(filenameselimination{i,1},annot_filename))
        a=0;
    end
end

end