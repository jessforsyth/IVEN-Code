function [cell_id,x,y,z,ch1,ch2,ch3,ch2_adj,ch3_adj]=import_2chan_data(table,headings);
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
    ch2_adj=table(headings:end,8);
    ch3_adj=table(headings:end,9);
    
end