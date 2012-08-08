function mouseLeftClick(objH,evt)
    global lineHandle returnedPath liveWireEngine imagePixels digitizedPath;
    keepGoing =0;
    if strcmp(get(objH,'SelectionType'),'alt') %Right click, stop digitizing
        keepGoing =0;
        if exist('returnedPath','var')
            digitizedPath = cat(1,digitizedPath,returnedPath+1);
        end
        disp('Stopped liveWire');
        set(gcf,'WindowButtonMotionFcn','');
        disp('Done digitizing');
        plot(digitizedPath(:,2),digitizedPath(:,1),'b-');
    else                    %Left click, set seedPoint
        if exist('returnedPath','var')
            digitizedPath = cat(1,digitizedPath,returnedPath+1);
        end
        %set(gcf,'WindowButtonMotionFcn','');    %Set windowButtonMotion off for the duration of calculations
        set(gcf,'WindowButtonMotionFcn',@mouseMoved);   %Turn on windowButtonMotionFcn
        point = get(get(objH,'Children'),'CurrentPoint');
        seedPoint = [round(point(1,1)), round(point(1,2))]; 
        %disp(['Seed ' num2str(seedPoint(:)')]);
        lineHandle = plot(seedPoint(1),seedPoint(2),'r-');
        if ~exist('digitizedPath','var')
            digitizedPath(1,1) = seedPoint(2);
            digitizedPath(1,2) = seedPoint(1);
        end
        drawnow;
        disp('Seed set');
        disp('Seed running');
        liveWireEngine.setSeed(seedPoint(2)-1,seedPoint(1)-1);    %N.B. row, column and index differs by 1 BTW matlab and JAVA!!!
        disp('Seed set ready');
    end
end