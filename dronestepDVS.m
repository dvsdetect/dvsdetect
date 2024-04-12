clear all; clc;

addpath(genpath('C:\Users\wenhao.lu\Desktop\code\30Aug_code_chip_testing\Matlab_AER_vision_functions'));
addpath(genpath('C:\Users\wenhao.lu\Desktop\code\30Aug_code_chip_testing\HIST_EXTRACTOR'))

wp = 55e3;%50e3

wo_entry_exit=1; % ??
extra_xSize=4; % ??
extra_ySize=3; % ??
max_x=346; % It's the size of the frame along x, max_x=86
max_y=260; % It's the size of the frame along y, max_y=64
max_hist_size=64; % max histogram size is taken 64
max_hist_amp = 32;
classes_to_augment = [2,5,6];
no_classes = 7;
n_objs = 8;
RP_type = "histogram_based_RP";
%RP_type = "CCA_based_RP";
com_path = 'C:/Users/wenhao.lu/Desktop/code/MINION_ANNOTATIONS';
Recording = 'All_TD_files2';
Annotation = 'Annotation';
file_type = '.mat';
% database_save_path = '/code/deepak/neuprowler_matlab/DATABASE';
% database_save_path = '/code/deepak/neuprowler_matlab/CCA_RP_DATABASE';
database_save_path = 'C:/Users/wenhao.lu/Desktop/code/MINION_ANNOTATIONS/HIST_RP_DATABASE';
dataset_prefixs = {...
%                   'sample1',...   %
                    'sample2',...          %
                   };
   
cs = Constants(com_path,Recording,Annotation,...
                     dataset_prefixs,file_type,wp,max_x,max_y,max_hist_size,...
                     no_classes, extra_xSize,extra_ySize,wo_entry_exit,...
                     max_hist_amp,classes_to_augment,...
                     database_save_path,n_objs);
                     
                 
                 
all_database = feature_database(cs);

for folder_id = cs.dataset_prefixs
    fid = cell2mat(folder_id);
    %ID = erase(annot_id,".txt");
    TD = cs.load_tracks_with_TD2(fid);
%     bin_image = all_database.explore_video_using_RPadd(TD,RP_type);
%     [record,current] = trackingalg3(bin_image);
%    [record,current] = all_database.explore_video_using_RPcombine(TD,RP_type); 
%    [record,current] = all_database.explore_video_using_RPCCLcombine(TD,RP_type);
%    [record,current] = all_database.explore_video_using_RPCCLcombineDVS(TD,RP_type);
    [record,current] = all_database.explore_video_using_RPCCLcombineDVSslens(TD,RP_type);  %%%%%wanted
%    [record,current] = all_database.explore_video_using_RPCCLcombineDVSslens2(TD,RP_type);
%    [record,current] = all_database.explore_video_using_RPCCLcombineDVSslensdownsample(TD,RP_type);
%    [record,current] = all_database.explore_video_using_RPCCLbased(TD,RP_type);  %%%%%wanted
%    [record,current] = all_database.explore_video_using_RPorg(TD,RP_type);  %%%%%wanted

end