function run_batch_IVEN
%Internal Versus External Neighbourhood (IVEN) quantification December 2020
%Jessica E. Forsyth- Plusa Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By calling this function, you run the batch version of 'IVEN'. No GUIs are
%used, and no manual checking of cell classification. 
%We strongly recommend that you analyse datasets using the original version
%of IVEN to ensure accurate classification of cells, however we provide
%this for more simplistic data sets that do not require such interrogation.
%

close all 
clear all
clc

%NOTE:
%-this code with automatically assign cells as inside or outside as based
%on the convex hull and shrink factor assigned below....
shrink=0; %if shrink =0 = basic convex hull, if shrink=1 = internalised CH

%Thresholding parameters, please change these as is appropriate to your
%data:
%Option=1 [Threshold using the inferred max nbr distance]
%Option=2 [Threshold set to value e.g. cell diameter]
%Option=3 [No thresholding of neighbours]
%Option=4 [Threshold using inferred max nbr distance for inside and outside
%cells] 
%.....if using Option 1 or 4 set the value of k
%params
option=4;
thresh_param=0.5;

%stype the number of headings rows your data has
headings=1;

[file,path]=uigetfile('*.*','MultiSelect','on');

if class(file)=='char'
    file=reshape(file,[size(file),1])
    file=cellstr(file);
end

num_files=length(file); 

for file_num=1:num_files
    fprintf('File %d out of %d : ',file_num,num_files)
    fprintf([file{file_num},'\n'])                 %this displays the filename for reference
    fname=fullfile(path,file{file_num});           %stitch file and path 
    table=table2array(readtable(fname));           %read in data from excel spreadsheet

    s=size(table);
    num_cells=s(1);                                 %calculate the number of cells in the dataset
    if num_cells<=4
        %If there are 1,2, 3 or 4 cells, there are too few cells to compute a
        %Voronoi diagram, it is also intuitive that these cells will have n-1
        %neighbours. Please do not run the code for these cases.
        fprintf('%s skipped, (too few cells to generate Voronoi diagram).\n',char(file))
        con=false;
    else 
        con=true;
    end

    %now extract data from imported excel files/
    if con==true
        cols=s(2);
        if cols==9
            'Green channel not present.';
            [cell_id,x,y,z,ch1,ch2,ch3,ch2_adj,ch3_adj]=import_2chan_data(table,headings);
            trms2=find(isnan(ch2_adj));
            ch2_adj(trms2)=0;
            trms3=find(isnan(ch3_adj));
            ch3_adj(trms3)=0;

            ch4=ones(1,length(x))';
            ch4_adj=ones(1,length(x))';
            fprintf(['Channel 4=false, only three channels present.','\n'])

        elseif cols==11
            %This part of the loop should always be used if data is input into the excel file correctly.
            % For ease we suggest filling 'empty' channels with 1's so it is
            % obvious that the channel data is not 'real' data, simply padding
            % data to allow for consistent input.
            'Green channel present. If this is not true, please check your file.';
            [cell_id,x,y,z,ch1,ch2,ch3,ch4,ch2_adj,ch3_adj,ch4_adj]=import_3chan_data(table,headings);
            trms2=find(isnan(ch2_adj));
            ch2_adj(trms2)=0;
            trms3=find(isnan(ch3_adj));
            ch3_adj(trms3)=0;
            trms4=find(isnan(ch4_adj));
            ch4_adj(trms4)=0;
        else
            error('Data not formatted correctly, incorrect number of columns, please reformat your data and try again.') 
        end
    end
    
    num_cells=length(x);
    fprintf('Number of cells - %d \n',num_cells)
    xyz=[x,y,z];
    TRI=delaunay(xyz);       %generate the delaunay triangulation (DT) of the imported points.
    
    num_tetra=length(TRI);   %calculate the number of tetrahedrons generated in the DT.

    [outside,inside]=outside_selection(cell_id,xyz,num_cells,shrink);  %now classify cells as inside or outside
    
    %now calculate the number of neighbours for the cells
    [nbrs,num_nbrs,nbrs_new,num_nbrs_new]=dt_nbr_calc(xyz,TRI,num_cells,num_tetra,option,thresh_param,outside); %calculate the number of 'true' neighbours

      
    %As an extra, and to demonstrate what sort of data IVEN can provide, we
    %added this part on, this describes the composition of the number of
    %neighbours that are outside cells. This was used to identify the mural
    %and polar TE cells. 
    %This part goes through the list of neighbours for each cell and sees how
    %many of them were outside. 
    %Similar analyses to this can be added into the code for output, either
    %into the excel file or into the command window. 
    nbr_comp=zeros(num_cells,4); %generates a table with cell ID | cell classification | #nbrs that are TE |#nbrs that are INS
    for c1=1:num_cells
        nbr_comp(c1,1)=cell_id(c1);
        nbr_comp(c1,2)=ismember(c1,outside);
        for c2=1:num_cells
            if nbrs_new(c1,c2)==1
                if ismember(c2,outside)==1
                    %count number of nbrs that are TE
                    nbr_comp(c1,3)=nbr_comp(c1,3)+1;
                else
                    %count number of nbrs that are ICM
                    nbr_comp(c1,4)=nbr_comp(c1,4)+1;
                end
            end
        end
    end
    
    %Another extra! The output of the cell_IDs of neighbours of each cell 
    [nbr_list]=nbr_list_gen(nbrs_new,cell_id,num_cells);
    
    %Final extra! (see how easy it is to add more modules to egt the most
    %out of the data. To output the median distance between a cell and its
    %neighbours. 
    [med_nbr_dist]=median_nbr_dist_calc(nbrs_new,num_cells,num_nbrs_new,xyz);   

    close all
    finalfig(cell_id,xyz,num_cells,ch2_adj,ch3_adj,ch4_adj,outside,inside,path,file{file_num}) %generate final figure of embryo (with cell IDs labelled and cell class. shown
    close all
    output_file_3chan(file{file_num},path,cell_id,x,y,z,ch1,ch2,ch3,ch4,ch2_adj,ch3_adj,ch4_adj,num_nbrs,num_nbrs_new,outside,num_cells,nbr_comp(:,3),med_nbr_dist,nbr_list); %create and save output file as .xls file.
    fprintf(['File saved.','\n'])
    fprintf('...................\n')
end
fprintf(['All files processed.','\n'])
    
end
    
        