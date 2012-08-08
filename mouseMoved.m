function mouseMoved(objH,evt)
    global lineHandle returnedPath liveWireEngine imagePixels digitizedPath;
    point = get(get(objH,'Children'),'CurrentPoint');
    targetPoint = [round(point(1,1)), round(point(1,2))];
    if targetPoint(1) < 1 || targetPoint(1) > size(imagePixels,2) || targetPoint(2) < 1 || targetPoint(2) > size(imagePixels,1)
        disp('Out of image');
    else
        disp(['Set target']);
%         whereFrom = liveWireEngine.whereFrom;
        returnedPath = liveWireEngine.returnPath(targetPoint(2)-1,targetPoint(1)-1); %N.B. row, column!!!
        
%         keyboard;
        if isempty(returnedPath)
            disp(['Not there yet']);
           return; 
        end
        disp(['Target set']);
        %Display image
        set(lineHandle,'XData',returnedPath(:,2)+1,'YData',returnedPath(:,1)+1); %N.B. row, column and index differs by 1!!!
        drawnow;
    end
end