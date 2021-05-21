function finalfig(cell_id,xyz,num_cells,ch1_int,ch2_int,ch3_int,outside,inside,path,file)
%This function plots the points of the embryo, displays their
%classification as inside or outside as before, and also maps the cell ID
%to a point. So tht through cross comparison with the output excel files,
%interesting features/anomalies can be identified within the embryo.

FigH = figure('position',[50 0 600 500],'visible','off');
movegui(FigH,'west')
axes('units','pixels', ...
     'position',[100 150 400 350], 'NextPlot', 'add');
axis('equal');

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);
markersize=120; 
allplh=scatter3(x,y,z,markersize,'o','MarkerFaceColor','b','MarkerEdgeColor','k');
hold on 

nums=num2str(cell_id);
nums2=cellstr(nums);
text(x,y,z,nums2,'Fontsize',10);
hold on

view(49,11)
grid on
hold on
shrink=0;
k=boundary(x,y,z,shrink);

plH = scatter3(x(outside),y(outside),z(outside),markersize,'o','MarkerFaceColor','w','MarkerEdgeColor','k');
hold on
trisurf(k,x,y,z,'FaceAlpha',0.05);

fname=fullfile(path,strcat('fig_',file,'.fig'));
savefig(gcf,char(fname))
    
end

