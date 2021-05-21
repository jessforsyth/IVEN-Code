function [med_nbr_dist]=median_nbr_dist_calc(nbrs_new,num_cells,num_nbrs_new,xyz)
%this function calculates the distance between a cell and all of its true
%neighbours and calculates the median value. 

med_nbr_dist=zeros(num_cells,1);
for cell=1:num_cells
    dist_temp=[];
    for cell2=1:num_cells
        if nbrs_new(cell,cell2)==1
            dist_temp=[dist_temp,pdist([xyz(cell,:);xyz(cell2,:)],'euclidean')]; %calculate the Euclidean distance between two neighbour
        end
    end
    med_nbr_dist(cell,1)=median(dist_temp);    
end

end
