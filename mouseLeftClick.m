function mouseLeftClick(objH,evt)
    global lineHandle parameters whereFrom keepGoing seedRunning;
    keepGoing =0;
    if strcmp(get(objH,'SelectionType'),'alt') %Right click, stop digitizing
        keepGoing =0;
        disp('Stopped liveWire');
        set(gcf,'WindowButtonMotionFcn','');
        disp('Done digitizing');
    else                    %Left click, set seedPoint
        %set(gcf,'WindowButtonMotionFcn','');    %Set windowButtonMotion off for the duration of calculations
        set(gcf,'WindowButtonMotionFcn',@mouseMoved);   %Turn on windowButtonMotionFcn
        point = get(get(objH,'Children'),'CurrentPoint');
        seedPoint = [round(point(1,2)), round(point(1,1))];
        %disp(['Seed ' num2str(seedPoint(:)')]);
        lineHandle = plot(seedPoint(2),seedPoint(1),'r*-');
        drawnow;
        parameters.seedPoint = seedPoint;
        disp('Seed set');
        while seedRunning == 1
           pause(0.1);
        end
        disp('Done waiting previous running seed');
        seedRunning = 1;
        keepGoing = 1;
        disp('Seed running');
        liveWireSetSeed;    %Set seed
        disp('Seed set ready');
        %disp(['Seed ' num2str(parameters.seedPoint(:)')]);
%         set(gcf,'WindowButtonMotionFcn',@mouseMoved);   %Turn on windowButtonMotionFcn
    end
end