classdef Constants 
   properties 
      window_period
      com_path
      rec_folder
      annot_folder
      dataset_prefixs
      extra_xSize
      extra_ySize
      maxX
      maxY
      max_hist_size
      file_type
      wo_entry_exit
      stride
      hist_consts
      RP_consts
      classes_to_augment = []
      database_save_path
      no_classes
   end
   methods
       function cs = Constants(com_path,rec_folder,annot_folder,...
                     dataset_prefixs,file_type,wp,maxX,maxY,max_hist_size,...
                     no_classes,extra_xSize,extra_ySize,wo_entry_exit,...
                     max_hist_amp,classes_to_augment,...
                     database_save_path,n_objs)
          if nargin > 0
             if ~exist(com_path,'dir')
                 error(['Folder not found: ',com_path]);
             else
                 cs.com_path  = com_path;
             end
             
             rec_f = [com_path,'/',rec_folder];
             if ~exist(rec_f,'dir')
                 error(['Folder not found: ',rec_f]);
             else
                 cs.rec_folder  = rec_f;
             end
             
             annot_f = [com_path,'/',annot_folder];
             if ~exist(rec_f,'dir')
                 error(['Folder not found: ',annot_f]);
             else
                 cs.annot_folder  = annot_f;
             end
             
             cs.dataset_prefixs = dataset_prefixs;
             cs.file_type = file_type;
             cs.window_period = wp;
             cs.maxX = maxX;
             cs.maxY = maxY;
             cs.no_classes = no_classes;
             cs.max_hist_size = max_hist_size;
             cs.extra_xSize = extra_xSize;
             cs.extra_ySize = extra_ySize;
             cs.wo_entry_exit = wo_entry_exit;
             cs.stride = 1;
             cs.hist_consts.th=1;
             cs.hist_consts.resolution=16;
             cs.hist_consts.max_hist_amp=max_hist_amp;
             cs.classes_to_augment = classes_to_augment;
             cs.database_save_path = database_save_path;

             if nargin>16
                cs.RP_consts.th = 1;
                cs.RP_consts.scale = 3;
                cs.RP_consts.n_objs = n_objs;
                cs.RP_consts.IoU = 0.1;
             end
          end
       end 
       
       
       function [TD, annot] = load_rec_file(cs,ID)
          TD = load([cs.rec_folder,'/TD_',ID,cs.file_type]);
          TD = TD.TD;
          annot = read_annotation([cs.annot_folder,'/',...
                  ID,'.txt']);
          annot.xSize = annot.xSize + cs.extra_xSize;
          annot.ySize = annot.ySize + cs.extra_ySize;
          
          time_offset=min(min(TD.ts),min(annot.ts))-10*cs.window_period;

          TD.ts=single(TD.ts-time_offset);
          annot.ts=single(annot.ts-time_offset);
       end

       function [TD,GT_annots] = load_tracks_with_TD(cs,ID)
          TD = load([cs.rec_folder,'/TD_',ID,cs.file_type]);
          TD = TD.TD;

          GT_annots = annotation.read_file([cs.annot_folder,'/',...
                  ID,'.txt']);

       end
   end
end