function [arrays]=intialise_arrays(num_files)

    arrays=[];
    %average number of neighbours for inside and outside cells per file
    arrays.avg_INS=zeros(num_files,1);
    arrays.avg_OUT=zeros(num_files,1);
    %number of neighbours of inside and outside cells from all files
    arrays.all_INS=[];
    arrays.all_OUT=[];
    %composition of neighbours for inside and outside cells
    arrays.comp_INS=[];  %percentage of neighbours that are also inside cells for inside cells
    arrays.comp_OUT=[];  %percentage of neighbours that are also outside cells for outside cells
    %numbers of neighbours of mural TE and polar TE (separately)
    arrays.mu_TE=[];
    arrays.po_TE=[];
    %numbers of cells per files
    arrays.cell_numbers=zeros(num_files,1);
    %number/percentage of mural TE cells per sample/embryo
    arrays.num_mural=[];
    arrays.perc_mural=[];
    
end