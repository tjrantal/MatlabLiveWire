function mouseMoved(objH,evt)
    global lineHandle returnedPath liveWireEngine imagePixels;
    point = get(get(objH,'Children'),'CurrentPoint');
    targetPoint = [round(point(1,1)), round(point(1,2))];
    if targetPoint(1) < 1 || targetPoint(1) > size(imagePixels,2) || targetPoint(2) < 1 || targetPoint(2) > size(imagePixels,1)
        disp('Out of image');
    else
        disp(['Set target']);
        returnedPath = liveWireEngine.returnPath(targetPoint(1),targetPoint(2));
        if isempty(returnedPath)
           return; 
        end
        disp(['Target set']);
        %Display image
        set(lineHandle,'XData',returnedPath(:,1),'YData',returnedPath(:,2));
        drawnow;
    end
end