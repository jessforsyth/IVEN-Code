function [outside,inside]=outside_selection(cell_id, xyz,num_cells,shrink)

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);
nums=num2str(cell_id);

%Now find outside cells (identified using Convex Hull) cells with a white outline, 
%and inside cells with a blue outline
k=boundary(x,y,z,shrink);  %convex hull
outside=[];
inside=[];
for i=1:num_cells
    if ismember(i,k)==1
        outside=[outside,i];
    else
        inside=[inside,i];
    end
end
            
end






