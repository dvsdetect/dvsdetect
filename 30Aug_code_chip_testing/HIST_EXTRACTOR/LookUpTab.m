classdef LookUpTab < handle
    properties 
        lookup 
        corner1x 
        corner1y
        corner2x
        corner2y
        
    end
    methods 
        function Table = LookUpTab(width,height)
            Table.lookup = zeros(floor(max(height,width)^2/4),1);    
            Table.corner1x = zeros(floor(max(height,width)^2/4),1);
            Table.corner1y = zeros(floor(max(height,width)^2/4),1);
            Table.corner2x = zeros(floor(max(height,width)^2/4),1);
            Table.corner2y = zeros(floor(max(height,width)^2/4),1);
        end
        function [out_lab,computes] = fillTab(Table,label,pixX,pixY)
            computes = 0;
            if (Table.lookup(label)==0)
                Table.lookup(label) = label;
                Table.corner1x(label) = pixX;
                Table.corner2x(label) = pixX;
                Table.corner1y(label) = pixY;
                Table.corner2y(label) = pixY;
                out_lab = label;
                computes = computes + 8;
            elseif (Table.lookup(label) == label)
                Table.corner1x(label) = min(pixX,Table.corner1x(label));
                Table.corner2x(label) = max(pixX,Table.corner2x(label));
                Table.corner1y(label) = min(pixY,Table.corner1y(label));
                Table.corner2y(label) = max(pixY,Table.corner2y(label));
                out_lab = label;
                computes = computes + 15;
            else
                while(Table.lookup(label) ~= label)
                    label = Table.lookup(label);
                    computes = computes + 3;
                end
                Table.corner1x(label) = min(pixX,Table.corner1x(label));
                Table.corner2x(label) = max(pixX,Table.corner2x(label));
                Table.corner1y(label) = min(pixY,Table.corner1y(label));
                Table.corner2y(label) = max(pixY,Table.corner2y(label));
                out_lab = label;
                computes =  computes +12;
            end
        end
        
        function [min_out,computes]= updateLookUp(Table,label1,label2)
            computes=0;
            if (label1<label2)
               while (Table.lookup(label1)~=label1)
                   label1 = Table.lookup(label1);
                   computes = computes + 3;
               end
               Table.lookup(label2) = label1;
               computes = computes + 1;
               if (Table.corner1x(label2)>0)
                   Table.corner1x(label1) = min(Table.corner1x(label1),...
                                                Table.corner1x(label2));
                   Table.corner2x(label1) = max(Table.corner2x(label1),...
                                                Table.corner2x(label2));
                   Table.corner1y(label1) = min(Table.corner1y(label1),...
                                                Table.corner1y(label2));
                   Table.corner2y(label1) = max(Table.corner2y(label1),...
                                                Table.corner2y(label2));
                   computes = computes + 10+8;
                   Table.emptylabel(label2);
                   computes = computes + 10;
               end
               min_out = label1;
            elseif (label2<label1)
%                 label2
               while (Table.lookup(label2)~=label2)
                   
                   label2 = Table.lookup(label2);
                   computes = computes + 3;
%                    label2
               end
               Table.lookup(label1) = label2;
               computes = computes + 1;
               if (Table.corner1x(label1)>0)
                   Table.corner1x(label2) = min(Table.corner1x(label1),...
                                                Table.corner1x(label2));
                   Table.corner2x(label2) = max(Table.corner2x(label1),...
                                                Table.corner2x(label2));
                   Table.corner1y(label2) = min(Table.corner1y(label1),...
                                                Table.corner1y(label2));
                   Table.corner2y(label2) = max(Table.corner2y(label1),...
                                                Table.corner2y(label2));
                   computes = computes + 10+8;
                   Table.emptylabel(label1);
                   computes = computes + 10;
               end
               min_out = label2;
            end
            computes = computes +1;
        end
        function emptylabel(Table,label)
            Table.corner1x(label) = 0;
            Table.corner2x(label) = 4*length(Table.corner1x) + 1;
            Table.corner1y(label) = 0 ;
            Table.corner2y(label) = 4*length(Table.corner1x) + 1;
        end
        
    end
end