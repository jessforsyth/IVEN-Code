function output_file_3chan(file,filepath,cell_id,x,y,z,ch1,ch2,ch3,ch4,ch2_adj,ch3_adj,ch4_adj,nbrs_num_old, nbrs_num_new,outside,num_cells,nbr_comp,med_nbr_dist,nbr_list)
%This function outputs data into an excel spreadsheet, separate to the
%input file. This way data cannot be overwritten, plus multiple comaprisons
%of data can be performed easily. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If extra columns to the output file need to be added:
%1. Ensure variable is input into the function (within the brackets after 'output_file_3chan(.....'
%2. Create a title in the list of column titles T, row 17
%3. Enter variable into T2, row 18

outside_id=zeros(num_cells,1);
for v=1:num_cells
    if ismember(v,outside)==1
        outside_id(v)=1;
    end
end


T=[string('Cell_ID'),string('X'),string('Y'),string('Z'),string('Ch1'),string('Ch2'),string('Ch3'),string('Ch4'),string('Ch2_adj'),string('Ch3_adj'),string('Ch4_adj'),string('Outside'),string('#nbrs_old'), string('#nbrs_new'), string('#nbrs_that_are_outside'),string('Median Dist to nbrs'),string('nbr_ID_list')];
T2=[cell_id,x,y,z,ch1,ch2,ch3,ch4,ch2_adj,ch3_adj,ch4_adj,outside_id,nbrs_num_old, nbrs_num_new,nbr_comp,med_nbr_dist,nbr_list];
T2=sortrows(T2,12);
T=[T;T2];

fname=fullfile(filepath,strcat('nbrs_',file));
[dir,name,ext]=fileparts(fname);
fname2=fullfile(dir,strcat(name,'.txt'));

%if using a windows machine, save data as an excel file
xlswrite(char(fname), T);

%if using a Mac, save data as a txt file and later open using Excel and delimiter
%{
fid=fopen(fname2,'wt')
for ii=1:size(T,1)
    if ii==1
        fprintf(fid,'%s ',T(ii,:));
        fprintf(fid,'\n');
    end
    
    fprintf(fid, '%.0f %f %f %f %f %f %f %f %f %f %f %.0f %.0f %.0f ', T(ii,:)); %need a valid format
    fprintf(fid, '\n');
end
%}

end