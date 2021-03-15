function [fig2,uic1,uid,uic2,uit,uic3,uic4,uid2]=thresholdChoice()

fig2=uifigure('Position', [360 398 560 220]);
drawnow

uil=uilabel(fig2,'Text','What threshold would you like to apply to this dataset?','Position', [20,200, 500,22],'Fontsize',13,'Fontweight','bold');

uic1=uicheckbox(fig2,'Text', 'Automatic (using the IQR) k=','Position', [20 170,250, 22],'Value',0);
uid=uidropdown(fig2,'Items',{'0.0','0.5','1.0','1.5','2.0','2.5','3.0'}, 'Value', '0.5','Position', [275, 170, 250, 22]);

uic2=uicheckbox(fig2,'Text', 'Pre-set value','Position', [20 140,250, 22]);
uit=uitextarea(fig2,'Value',num2str(30.000),'Position', [275 140,250, 22]);

uic3=uicheckbox(fig2,'Text','No threshold', 'Position', [20 110,250, 22]);

uic4=uicheckbox(fig2,'Text', 'Outside vs Inside automatic method k=','Position', [20 80,250, 22],'Value',1);
uid2=uidropdown(fig2,'Items',{'0.0','0.5','1.0','1.5','2.0','2.5','3.0'}, 'Value', '0.5','Position', [275, 80, 250, 22]);


uib=uibutton(fig2,'state','Text', 'OK', 'Position', [480 20 50 22]);

uic1.ValueChangedFcn=@(uic1,event) option1(event,fig2,uic1,uic2,uic3,uic4);
uic2.ValueChangedFcn=@(uic2,event) option2(event,fig2,uic1,uic2,uic3,uic4);
uic3.ValueChangedFcn=@(uic3,event) option3(event,fig2,uic1,uic2,uic3,uic4);
uic4.ValueChangedFcn=@(uic4,event) option4(event,fig2,uic1,uic2,uic3,uic4);


waitfor(uib,'Value');

    function option1(event,fig2,uic1,uic2,uic3,uic4);
        drawnow
        if event.Value==1
            uic2.Value=0;
            uic3.Value=0;
            uic4.Value=0;
        end
    end

    function option2(event,fig2,uic1,uic2,uic3,uic4);
        drawnow
        if event.Value==1
            uic1.Value=0;
            uic3.Value=0;
            uic4.Value=0;
        end
    end

    function option3(event,fig2,uic1,uic2,uic3,uic4);
        drawnow
        if event.Value==1
            uic1.Value=0;
            uic2.Value=0;
            uic4.Value=0;
        end
    end

    function option4(event,fig2,uic1,uic2,uic3,uic4);
        drawnow
        if event.Value==1
            uic1.Value=0;
            uic2.Value=0;
            uic3.Value=0;
        end
    end


end
