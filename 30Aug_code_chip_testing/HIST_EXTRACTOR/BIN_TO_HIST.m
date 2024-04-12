clear all; clc;
 
addpath(genpath(['/home/deepak/STEE-Tracker/',...
    'NTU_MATLAB_tracker_code/Matlab_AER_vision_functions']));

wp = 66e3;
wo_entry_exit=1; % ??
extra_xSize=4; % ??
extra_ySize=3; % ??
max_x=240; % It's the size of the frame along x
max_y=180; % It's the size of the frame along y
max_hist_size=64; % max histogram size is taken 64
max_hist_amp = 32;
classes_to_augment = [2,5,6];
no_classes = 7;
n_objs = 16;
RP_type = "histogram_based_RP";
% RP_type = "CCA_based_RP";
com_path = '/media/NEUPROWLER/common/MINION_ANNOTATIONS';
Recording = 'All_TD_files';
Annotation = 'Annotation';
file_type = '.mat';
% database_save_path = '/code/deepak/neuprowler_matlab/DATABASE';
% database_save_path = '/code/deepak/neuprowler_matlab/CCA_RP_DATABASE';
database_save_path = '/code/deepak/neuprowler_matlab/HIST_RP_DATABASE';
dataset_prefixs = {...
                    '20180711_ENG_3pm_12mm',...
                           '20180712_LT4_3pm_6mm_04',...
                           '20180713_ENG_4pm_12mm',...
                           '20180717_LT4_4pm_6mm',...
%                            '20180720_YIH_6pm_45mm',...
%                            '20180723_ENG_3pm_12mm',...
%                            '20180723_ENG_3pm_16mm',...
%                            '20180723_ENG_4pm_6_1mm',...
%                            '20180723_ENG_4pm_6mm',...
                   };
   
cs = Constants(com_path,Recording,Annotation,...
                     dataset_prefixs,file_type,wp,max_x,max_y,max_hist_size,...
                     no_classes, extra_xSize,extra_ySize,wo_entry_exit,...
                     max_hist_amp,classes_to_augment,... 
                     database_save_path,n_objs);
%% Start of Script



%% For each Folder [Train Test Validate]
% Run 1
all_database = feature_database(cs);
annot_filenames = dir([cs.annot_folder,'/','*.txt']);
annot_filenames = {annot_filenames.name};
for folder_id = cs.dataset_prefixs
    fid = cell2mat(folder_id);
    for annot_file = annot_filenames
        annot_id = cell2mat(annot_file);
        if startsWith(annot_id,fid)
            ID = erase(annot_id,".txt");
            [TD, annot_GT] = cs.load_rec_file(ID);
            all_database.explore_video_usign_GT(TD,annot_GT);
            save_filename = all_database.save_database(ID);
            disp(['Database saved at: ', save_filename]);
        end
        if (all_database.files_explored==15)
            break;
        end
    end
    if (all_database.files_explored==15)
        break;
    end
end


%% Variables while loading data again
del_class = 7;
human_class = 4;
classes_to_merge = [1,3];

%% Loading all mat train files

all_database = feature_database(cs);
train_folder = [cs.database_save_path,'/train_data'];
train_file_names = [train_folder,'/*.mat'];
feature_files =  dir(train_file_names);
feature_files = {feature_files.name};
for feat_id  = feature_files
    fid = cell2mat(feat_id);
    load([train_folder,'/',fid]);
    all_database.load_file(arr_to_save);
end

%% Updating the database
all_database.remove_class(del_class);
%%
all_database.remove_class(human_class);
%%
all_database.merge_classes(classes_to_merge);

%%
all_database.balance_DB(10000);

%% Saving train data in mat 
save_filename = all_database.save_database('train_file');
disp(['Database saved at: ', save_filename]);


%% Loading all mat test files
all_database = feature_database(cs);
train_folder = [cs.database_save_path,'/test_data'];
train_file_names = [train_folder,'/*.mat'];
feature_files =  dir(train_file_names);
feature_files = {feature_files.name};
for feat_id  = feature_files
    fid = cell2mat(feat_id);
    load([train_folder,'/',fid]);
    all_database.load_file(arr_to_save);
end

%% Removing the augmented data
all_database.remove_augmented_dat();
%% Updating the database
all_database.remove_class(del_class);
%%
all_database.remove_class(human_class);
%%
all_database.merge_classes(classes_to_merge);

%% Saving test data in mat
%%
save_filename = all_database.save_database('test_file');
disp(['Database saved at: ', save_filename]);








%% Making DB from RP

all_database = feature_database(cs);
annot_filenames = dir([cs.annot_folder,'/','*.txt']);
annot_filenames = {annot_filenames.name};
for folder_id = cs.dataset_prefixs
    fid = cell2mat(folder_id);
    for annot_file = annot_filenames
        annot_id = cell2mat(annot_file);
        if startsWith(annot_id,fid)
            ID = erase(annot_id,".txt");
            [TD, annot_GT] = cs.load_tracks_with_TD(ID);
%             all_database.save_comparison_video(TD,annot_GT);
            all_database.explore_video_using_RP(TD,annot_GT,...
                                            RP_type);
            save_filename = all_database.save_database(ID);
            disp(['Database saved at: ', save_filename]);
        end
        if (all_database.files_explored==15)
            break;
        end
    end
    if (all_database.files_explored==15)
        break;
    end
end

%% Variables while loading data again
del_class = 7;
human_class = 4;
classes_to_merge = [1,3];

%% Loading all mat train files

all_database = feature_database(cs);
train_folder = [cs.database_save_path,'/train_data'];
train_file_names = [train_folder,'/*.mat'];
feature_files =  dir(train_file_names);
feature_files = {feature_files.name};
for feat_id  = feature_files
    fid = cell2mat(feat_id);
    load([train_folder,'/',fid]);
    all_database.load_file(arr_to_save);
end

%% Updating the database
all_database.remove_class(del_class);
%%
all_database.remove_class(human_class);
%%
all_database.merge_classes(classes_to_merge);

%%
all_database.balance_DB(10000);

%% Saving train data in mat 
save_filename = all_database.save_database('train_file');
disp(['Database saved at: ', save_filename]);

%% Loading all mat test files
all_database = feature_database(cs);
train_folder = [cs.database_save_path,'/test_data'];
train_file_names = [train_folder,'/*.mat'];
feature_files =  dir(train_file_names);
feature_files = {feature_files.name};
for feat_id  = feature_files
    fid = cell2mat(feat_id);
    load([train_folder,'/',fid]);
    all_database.load_file(arr_to_save);
end

%% Removing the augmented data
all_database.remove_augmented_dat();
%% Updating the database
all_database.remove_class(del_class);
%%
all_database.remove_class(human_class);
%%
all_database.merge_classes(classes_to_merge);

%% Saving test data in mat
%%
save_filename = all_database.save_database('test_file');
disp(['Database saved at: ', save_filename]);