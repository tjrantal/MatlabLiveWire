function mouseMoved(objH,evt)
    global lineHandle parameters imagePixels whereFrom;
    point = get(get(objH,'Children'),'CurrentPoint');
    targetPoint = [round(point(1,2)), round(point(1,1))];
    if targetPoint(1) < 1 || targetPoint(1) > size(imagePixels,1) || targetPoint(2) < 1 || targetPoint(2) > size(imagePixels,2)
        disp('Out of image');
    else
        disp(['Set target']);
        ok = liveWireSetTarget(targetPoint);
        if ok == 0
           return; 
        end
        disp(['Target set']);
        %Display image
        set(lineHandle,'XData',parameters.returnedPath(2,:),'YData',parameters.returnedPath(1,:));
        drawnow;
    end
end