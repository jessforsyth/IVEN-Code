function [nbr_list]=nbr_list_gen(nbrs_new,cell_id,num_cells)

%this function will identify the cell IDs of each cell's neighbours for
%output into the nbrs_output file. 
nbr_list=strings(num_cells, 1) ;
for cell=1:num_cells
    nbrs_temp=find(nbrs_new(cell,:)==1); %find out which cells are neighbours. 
    cell_id_list=cell_id(nbrs_temp); %find the corresponding cell IDs
    nbrs_string=sprintf('%d, ', cell_id_list); % create string list delimited by commas
    nbr_list(cell,1)=nbrs_string(1:end-1);
end

    