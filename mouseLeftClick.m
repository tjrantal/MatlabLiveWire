function mouseLeftClick(objH,evt)
    global lineHandle returnedPath liveWireEngine imagePixels;
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
        seedPoint = [round(point(1,1)), round(point(1,2))]; 
%         keyboard;
        %disp(['Seed ' num2str(seedPoint(:)')]);
        lineHandle = plot(seedPoint(1),seedPoint(2),'r*-');
        drawnow;
        disp('Seed set');
        disp('Seed running');
        liveWireEngine.setSeed(seedPoint(2),seedPoint(1));    %N.B. row, column!!!
        disp('Seed set ready');
    end
end