function [nbrs,num_nbrs,nbrs_new,num_nbrs_new]=dt_nbr_calc(xyz,TRI,num_cells,num_tetra,option,thresh_param,outside)
%This function now caluclates and corrects the number of neighbours for
%each cell using the DT and cell diameter measurements. This function
%currently uses cell diameters relevant to the preimplantation mouse embryo
%but this can be changed for whichever system is being analysed. 


nbrs=zeros(num_cells,num_cells);
%using the previously caluclated DT, we now identify which cells are
%neighbours of each other and store this information in 'nbrs' 
for cell=1:num_cells
    for row=1:num_tetra
        
        if find(TRI(row,:)==cell)~=0
            for j=1:4 %(4 points to a tetrahedron)
                if j~=find(TRI(row,:)==cell)
                    nbrs(cell,TRI(row,j))=1;
                end
            end
        end
    end
end

num_nbrs=sum(nbrs,2); %count how many neighbours each cell has in total

%Now we correct the calculation of neighbours and account for the potrntail
%presence of a cavity.

%New method of identifying threshold of neighbour distance

if option==1
    tot_nbrs=sum(num_nbrs);

    D=zeros(1,tot_nbrs);
    nbrnum=0;
    for cell1=1:num_cells   %nested for loop, if there is a large number of cells, this may become slow- keep an eye out!
        for cell2=1:num_cells
            if nbrs(cell1,cell2)==1
                nbrnum=nbrnum+1;
                d=pdist([xyz(cell1,:);xyz(cell2,:)],'euclidean'); %calculate the Euclidean distance between two neighbour
                D(1,nbrnum)=d;
            end
        end
    end

    p75=prctile(D,75);           %calulate the 75th percentile
    iqrange=iqr(D);              %calculate the Inter Quartile Range
    k=thresh_param;                       %this can be changed if needed!
    threshold=p75+(k*iqrange);    %Calculate threshold
    
    threshold_list=ones(1,num_cells)*threshold;
    
elseif option==2
    threshold=thresh_param;
    threshold_list=ones(1,num_cells)*threshold;
    
elseif option==3
    tot_nbrs=sum(num_nbrs);

    D=zeros(1,tot_nbrs);
    nbrnum=0;
    for cell1=1:num_cells   %nested for loop, if there is a large number of cells, this may become slow- keep an eye out!
        for cell2=1:num_cells
            if nbrs(cell1,cell2)==1
                nbrnum=nbrnum+1;
                d=pdist([xyz(cell1,:);xyz(cell2,:)],'euclidean'); %calculate the Euclidean distance between two neighbour
                D(1,nbrnum)=d;
            end
        end
    end
    
    threshold=max(D)*100;
    threshold_list=ones(1,num_cells)*threshold;
    
elseif option==4
    tot_nbrs=sum(num_nbrs);
    
    D_te=zeros(1,tot_nbrs);
    D_icm=zeros(1,tot_nbrs);
    
    tenbr=0;
    icmnbr=0;
    
    nbrnum=0;
    for cell1=1:num_cells   %nested for loop, if there is a large number of cells, this may become slow- keep an eye out!
        for cell2=1:num_cells
            if nbrs(cell1,cell2)==1
                nbrnum=nbrnum+1;
                d=pdist([xyz(cell1,:);xyz(cell2,:)],'euclidean'); %calculate the Euclidean distance between two neighbour
                if ismember(cell1,outside)==1
                    tenbr=tenbr+1;
                    D_te(1,tenbr)=d;
                else
                    icmnbr=icmnbr+1;
                    D_icm(1,icmnbr)=d;
                end
            end
        end
    end
    
    D_te=D_te(1:tenbr);
    D_icm=D_icm(1:icmnbr);
    
    p75te=prctile(D_te,75);           %calulate the 75th percentile
    iqrangete=iqr(D_te);              %calculate the Inter Quartile Range
    k=thresh_param;                       %this can be changed if needed!
    thresholdte=p75te+(k*iqrangete);    %Calculate threshold
    
    p75icm=prctile(D_icm,75);           %calulate the 75th percentile
    iqrangeicm=iqr(D_icm);              %calculate the Inter Quartile Range
    k=thresh_param;                       %this can be changed if needed!
    thresholdicm=p75icm+(k*iqrangeicm);    %Calculate threshold
    
    %{
    threshold_list=ones(1,num_cells);
    for i=1:num_cells
        if ismember(i,TE)==1
            threshold_list(1,i)=thresholdte;
        else
            threshold_list(1,i)=thresholdicm;
        end
    end
    %}
    
end

if option==1 || option==2 || option==3
    fprintf(['Threshold-',num2str(threshold),'\n'])
else
    fprintf(['Threshold TE-',num2str(thresholdte),'\n'])
    fprintf(['Threshold ICM-',num2str(thresholdicm),'\n'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%here we now check the distances between neighbours and remove those greater than our threshold
nbrs_new=nbrs;

for cell1=1:num_cells   %nested for loop, if there is a large number of cells, this may become slow- keep an eye out!
    for cell2=1:num_cells
        if nbrs(cell1,cell2)==1
            d=pdist([xyz(cell1,:);xyz(cell2,:)],'euclidean'); %calculate the Euclidean distance between two neighbours
            
            if option==1 || option==2 ||option==3
                if d>threshold %compare to threshold for cell1 (i.e. if cell1=icm, use icm threshold)
                    nbrs_new(cell1,cell2)=0; %if d>threshold remove cell as neighbour
                end
            else
                if ismember(cell1,outside)==0 || ismember(cell2,outside)==0  %if one of the cells being compared is ICM
                    if d>thresholdicm                              %use ICM threshold
                        nbrs_new(cell1,cell2)=0;
                    end
                else
                    if d>thresholdte
                        nbrs_new(cell1,cell2)=0;
                    end
                end
            end
        end
    end
end

num_nbrs_new=sum(nbrs_new,2); %updated and corrected numbers of neighbours

            
end
                    