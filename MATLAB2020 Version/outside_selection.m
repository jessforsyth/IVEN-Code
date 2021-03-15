function [shrink,fig,uis,outside_original]=outside_selection(cell_id, xyz,num_cells,ch1_int,ch2_int,ch3_int)

fig=uifigure('Position',get(0,'Screensize'));
drawnow
fig.WindowState= 'max';

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);
nums=num2str(cell_id);

%Now find outside cells (identified using Convex Hull) cells with a white outline, 
%and inside cells with a blue outline
shrink=0;
k=boundary(x,y,z,shrink);  %convex hull
outside=[];
inside=[];
colour_list=ones(num_cells,3).*[0,0,1];
for i=1:num_cells
    if ismember(i,k)==1
        outside=[outside,i];
        colour_list(i,:)=[1,1,1];
    else
        inside=[inside,i];
    end
end
setappdata(fig,'Outside',outside);
setappdata(fig,'Inside',inside);

outside_original=outside;
inside_original=inside;

markersize=120; %must be even
%Plot cell centres/nuclei as a 3D scatter plot
uia=uiaxes(fig,'Position',[50 50 600 600],'XGrid',1,'YGrid',1,'ZGrid',1);
plt=scatter3(uia,x,y,z,markersize,'o','Linewidth',2);%,'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.3);
hold(uia,'on');
%text(uia,x,y,z,cellstr(nums),'Fontsize',11,'Color','r');
%hold(uia,'on');
trisurf(k,x,y,z,'FaceAlpha',0.0,'Parent',uia);
hold(uia,'on')
plt.CData=colour_list;



%generate colour schemes
ch1_int_n=(100*(round(ch1_int./max(ch1_int),2)));
C1=summer(101);
green=C1(abs(int16(ch1_int_n)-100)+1,:);
pltgreen=scatter3(uia,x,y,z,markersize/2,green,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',1.0);
set(pltgreen,'Visible','on')
hold(uia,'on');

ch2_int_n=(100*(round(ch2_int./max(ch2_int),2)));
C2=autumn(101); %red channel- %red=high, yellow=low
red=C2(abs(int16(ch2_int_n)-100)+1,:);
pltred=scatter3(uia,x,y,z,markersize/2,red,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',1.0);
set(pltred,'Visible','off')
hold(uia,'on');

ch3_int_n=(100*(round(ch3_int./max(ch3_int),2)))+1;
C3=gray(101); %far red channel- white=high, black=low
fred=C3(int16(ch3_int_n),:);
pltfred=scatter3(uia,x,y,z,markersize/2,fred,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',1.0);
set(pltfred,'Visible','off')
hold(uia,'on');


setappdata(fig,'colour_list',colour_list);
colour_list_original=colour_list;

%button to start manual correction of outside cells
uib=uibutton(fig, 'push', 'Text', 'Manual Cell Classification', 'Position', [700, 500, 250, 22],'Fontsize',13);
uib2=uibutton(fig, 'push', 'Text', 'Reset', 'Position', [700, 425, 250, 22],'Fontsize',13);
uib3=uibutton(fig, 'push', 'Text', 'Outside', 'Position', [975, 500, 120, 22],'Fontsize',13);
uib4=uibutton(fig, 'push', 'Text', 'Inside', 'Position', [1100, 500, 120, 22],'Fontsize',13);

uib5=uibutton(fig,'push','Text', 'Shrink Factor- 0.0000', 'Position',[700 300 350 22],'fontsize',13);
uis=uislider(fig,'Position', [700 275 350 3],'Value', 0, 'Limits', [0,1.0],'MinorTicks',[]);

uib6=uibutton(fig,'push','Text', 'Rotate', 'Position',[975 425 250 22],'fontsize',13);

uic1=uicheckbox(fig,'Text', 'Ch2', 'Position', [700, 150,50,22],'Fontsize',13,'Value',1);
uic2=uicheckbox(fig,'Text', 'Ch3', 'Position', [700, 120,50,22],'Fontsize',13);
uic3=uicheckbox(fig,'Text', 'Ch4', 'Position', [700, 90,50,22],'Fontsize',13);

uib7=uibutton(fig,'state','Text','Finish', 'Position', [975, 90,200,70],'Fontsize',30);

uib.ButtonPushedFcn=@(uib,event) manualCorrection(fig, uia);
uib2.ButtonPushedFcn=@(uib2,evet) resetCorrection(fig, uia, uis, plt, outside_original,inside_original,colour_list_original);
uib3.ButtonPushedFcn=@(uib3,event) addOutside(fig,uia, plt,uic1,uic2,uic3);
uib4.ButtonPushedFcn=@(uib4,event) addInside(fig,uia, plt,uic1,uic2,uic3);

uis.ValueChangingFcn=@(uis,event) shrinkChange(event,fig, uia, uib5, plt,x,y,z,num_cells);

uib6.ButtonPushedFcn=@(uib6,event) rotateButton(event,fig,uia);

uic1.ValueChangedFcn=@(uic1,event) channel1Show(event,fig,pltgreen,pltred,pltfred,uic1,uic2,uic3);
uic2.ValueChangedFcn=@(uic2,event) channel2Show(event,fig,pltgreen,pltred,pltfred,uic1,uic2,uic3);
uic3.ValueChangedFcn=@(uic3,event) channel3Show(event,fig,pltgreen,pltred,pltfred,uic1,uic2,uic3);

%uib7.ButtonPushedFcn=@(uib7,event) finishCorrection(fig, uia);

waitfor(uib7,'Value');

    %function finishCorrection(fig,uia)  
     %   pause('off')
    %end
        

    function manualCorrection(fig, uia)
        drawnow
        brush(uia,'on')
    end

    function resetCorrection(fig, uia, uis, plt, outside_original,inside_original,colour_list_original)
        setappdata(fig,'Outside',outside_original);
        setappdata(fig,'Inside',inside_original);
        setappdata(fig,'colour_list',colour_list_original);
        plt.CData=colour_list_original;
        uis.Value=0;
    end

        

    function addOutside(fig, uia, plt,uic1,uic2,uic3)
        drawnow
        colour_list=getappdata(fig,'colour_list');
        outside=getappdata(fig,'Outside');
        inside=getappdata(fig,'Inside');
        if isempty(find(get(plt,'BrushData')==1))==0  || isempty(find(get(pltgreen,'BrushData')==1))==0 || isempty(find(get(pltred,'BrushData')==1))==0 || isempty(find(get(pltfred,'BrushData')==1))==0
            if isempty(find(get(plt,'BrushData')==1))==0 
                brush_val=find(get(plt,'BrushData')==1);
            elseif uic1.Value==1 && isempty(find(get(pltgreen,'BrushData')==1))==0;
                brush_val=find(get(pltgreen,'BrushData')==1);
            elseif uic2.Value==1 && isempty(find(get(pltred,'BrushData')==1))==0;
                brush_val=find(get(pltred,'BrushData')==1);
            elseif uic3.Value==1 && isempty(find(get(pltfred,'BrushData')==1))==0;
                brush_val=find(get(pltfred,'BrushData')==1);
            end
            
            outside=[outside,brush_val];
            inside(find(inside==brush_val))=[];
            colour_list(brush_val,:)=[1 1 1];          
            plt.CData=colour_list;
            setappdata(fig,'colour_list',colour_list);
            setappdata(fig,'Outside',outside);
            setappdata(fig,'Inside',inside);
        end
    end

    function addInside(fig, uia, plt,uic1,uic2,uic3)
        drawnow
        colour_list=getappdata(fig,'colour_list');
        outside=getappdata(fig,'Outside');
        inside=getappdata(fig,'Inside');
        if isempty(find(get(plt,'BrushData')==1))==0  || isempty(find(get(pltgreen,'BrushData')==1))==0 || isempty(find(get(pltred,'BrushData')==1))==0 || isempty(find(get(pltfred,'BrushData')==1))==0
            if isempty(find(get(plt,'BrushData')==1))==0 
                brush_val=find(get(plt,'BrushData')==1);
            elseif uic1.Value==1 && isempty(find(get(pltgreen,'BrushData')==1))==0;
                brush_val=find(get(pltgreen,'BrushData')==1);
            elseif uic2.Value==1 && isempty(find(get(pltred,'BrushData')==1))==0;
                brush_val=find(get(pltred,'BrushData')==1);
            elseif uic3.Value==1 && isempty(find(get(pltfred,'BrushData')==1))==0;
                brush_val=find(get(pltfred,'BrushData')==1);
            end
            
            inside=[inside,brush_val];
            outside(find(outside==brush_val))=[];
            colour_list(brush_val,:)=[0 0 1];          
            plt.CData=colour_list;
            setappdata(fig,'colour_list',colour_list);
            setappdata(fig,'Outside',outside);
            setappdata(fig,'Inside',inside);
        end
    end

    function shrinkChange(event,fig, uia, uib5, plt,x,y,z,num_cells)
        drawnow
        
        shrink=event.Value;
        k=boundary(x,y,z,shrink);
        outside=[];
        inside=[];
        
        colour_list=ones(num_cells,3).*[0,0,1];
        for i=1:num_cells
            if ismember(i,k)==1
                outside=[outside,i];
                colour_list(i,:)=[1,1,1];
            else
                inside=[inside,i];
            end
            
        end
        
        setappdata(fig,'colour_list',colour_list);
        setappdata(fig,'Outside',outside);
        setappdata(fig,'Inside',inside);  
        
        plt.CData=colour_list;
        txt=sprintf('%s %.4f','Shrink Factor-',round(event.Value,4));
        set(uib5,'Text',txt);
    end

    function rotateButton(event,fig,uia)
        rotate3d(fig,'on')        
    end

        
        
    function channel1Show(event,fig,pltgreen,pltred,pltfred,uic1,uic2,uic3)
        if event.Value==1
            set(pltgreen,'Visible','on');
            set(pltred,'Visible','off');
            set(pltfred,'Visible','off');
            uic2.Value=0;
            uic3.Value=0;
        else
            set(pltgreen,'Visible','off');
        end
    end

    function channel2Show(event,fig,pltgreen,pltred,pltfred,uic1,uic2,uic3)
        if event.Value==1
            set(pltgreen,'Visible','off');
            set(pltred,'Visible','on');
            set(pltfred,'Visible','off');           
            uic1.Value=0;
            uic3.Value=0;
        else
            set(pltred,'Visible','off');
        end
    end

    function channel3Show(event,fig,pltgreen,pltred,pltfred,uic1,uic2,uic3)
        if event.Value==1
            set(pltgreen,'Visible','off');
            set(pltred,'Visible','off');
            set(pltfred,'Visible','on');            
            uic1.Value=0;
            uic2.Value=0;
        else
            set(pltfred,'Visible','off');
        end
    end

            

end






