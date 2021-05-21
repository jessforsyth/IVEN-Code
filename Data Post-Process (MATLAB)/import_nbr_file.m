function [cell_id,x,y,z,ch1,ch2,ch3,ch4,ch2_adj,ch3_adj,ch4_adj,outside,nbrs_old,nbrs_new,nbr_comp,med_nbr_dist]=import_nbr_file_python(table,headings);
% this function extracts your data into the various variables. 
%function str2double simply converts the values in the table into 'usable'
%data. 
    cell_id=table(headings:end,1);
    x=table(headings:end,2);
    y=table(headings:end,3);
    z=table(headings:end,4);
    ch1=table(headings:end,5);
    ch2=table(headings:end,6);
    ch3=table(headings:end,7);
    ch4=table(headings:end,8);
    ch2_adj=table(headings:end,9);
    ch3_adj=table(headings:end,10);
    ch4_adj=table(headings:end,11);
    outside=table(headings:end,12);
    nbrs_old=table(headings:end,13);
    nbrs_new=table(headings:end,14);
    nbr_comp=table(headings:end,15);
    med_nbr_dist= table(headings:end,16);
    
    
end