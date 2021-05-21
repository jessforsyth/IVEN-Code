close all 
clear all
clc

%This script generates embryos schematics with mural and polar TE cells
%colour coded. 
%type the number of headings rows your data has
headings=1;

[file,path]=uigetfile('*.*','MultiSelect','on'); %This line opens a file directory window to allow you to select the files you would like to copile/post-process

if class(file)=='char'  %if you have selected multiple files, then file with be in format char, must reshape
    file=reshape(file,[size(file),1]);
    file=cellstr(file);
end
num_files=length(file); 

[arrays]=intialise_arrays(num_files);

warningid='MATLAB:table:ModifiedAndSavedVarnames';
warning('off',warningid);

for file_num=1:num_files
    fprintf('File %d out of %d : ',file_num,num_files)
    fprintf([file{file_num},'\n'])                 %this displays the filename for reference
    fname=fullfile(path,file{file_num});           %stitch file and path 
    table=table2array(readtable(fname),'PreserveVariableNames','True'); %read in data from excel spreadsheet
    s=size(table);
    num_cells=s(1);                      %calculate the number of cells in the dataset
    
    arrays.cell_numbers(file_num,1)=num_cells;
    
    [~,order]=sort(table(:,1));
    table=table(order,:);
    
    %Python extraction  
    %if you are extracting data from files generated by python use this
    %[~,~,~,~,~,~,~,~,~,~,~,OUT,nbrs_new,nbr_comp]=import_nbr_file_python(table,headings);
    
    %matlab extraction
    %if you are extracting data from files generated by matlab use this
    [~,x,y,z,~,~,~,~,~,~,~,OUT,nbrs_old,nbrs_new,nbr_comp]=import_nbr_file(table,headings);
    
    ins_nbrs=[];
    out_nbrs=[];
    num_mu=0;
    num_out=0;
    
    x_ins=[];
    y_ins=[];
    z_ins=[];
    x_mu=[];
    y_mu=[];
    z_mu=[];
    x_po=[];
    y_po=[];
    z_po=[];
    
    for i=1:num_cells
        if OUT(i)==0
            arrays.all_INS=[arrays.all_INS;nbrs_new(i)];
            ins_nbrs=[ins_nbrs;nbrs_new(i)];
            arrays.comp_INS=[arrays.comp_INS;(100*((nbrs_new(i)-nbr_comp(i))/nbrs_new(i)))];
            x_ins=[x_ins;x(i)];
            y_ins=[y_ins;y(i)];
            z_ins=[z_ins;z(i)];
        else
            nuout=num_out+1;
            arrays.all_OUT=[arrays.all_OUT;nbrs_new(i)];
            out_nbrs=[out_nbrs;nbrs_new(i)];
            arrays.comp_OUT=[arrays.comp_OUT;(100*(nbr_comp(i)/nbrs_new(i)))];

            if (nbr_comp(i)/nbrs_new(i))==1
                num_mu=num_mu+1;
                arrays.mu_TE=[arrays.mu_TE;nbrs_new(i)];
                x_mu=[x_mu;x(i)];
                y_mu=[y_mu;y(i)];
                z_mu=[z_mu;z(i)];
            else
                arrays.po_TE=[arrays.po_TE;nbrs_new(i)];
                x_po=[x_po;x(i)];
                y_po=[y_po;y(i)];
                z_po=[z_po;z(i)];
            end

        end
    end
 
    
    arrays.num_mural=[arrays.num_mural;num_mu];
    arrays.perc_mural=[arrays.perc_mural;100*(num_mu/num_out)];

    arrays.avg_INS(file_num,1)=mean(ins_nbrs);
    arrays.avg_OUT(file_num,1)=mean(out_nbrs);
   
    %generate figure showing the mural and polar TE cells   
    k=boundary([x_mu;x_po],[y_mu;y_po],[z_mu;z_po],0);
    
    
    figure()
    scatter3(x_ins,y_ins,z_ins,'ob','filled');
    hold on 
    scatter3(x_mu,y_mu,z_mu,'om','filled','MarkerEdgeColor','k');
    hold on
    scatter3(x_po,y_po,z_po,'oc','filled','MarkerEdgeColor','k');
    hold on
    trisurf(k,[x_mu;x_po],[y_mu;y_po],[z_mu;z_po],'FaceAlpha',0.0)
    title(file,'Interpreter','none')
    legend('ICM','Mural TE','Polar TE');
    axis equal
   
end

