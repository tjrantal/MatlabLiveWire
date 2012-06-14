function testJava
    
    close all;
    clear all;
    javaaddpath('.');
    global lineHandle returnedPath liveWireEngine imagePixels;
    dicomFileIn = 'IM000~10';
    imageInfo = dicominfo(dicomFileIn);
    imagePixels = double(dicomread(imageInfo));
    liveWireEngine = javaEngineLiveWire.LiveWireCosts(reshape(imagePixels',1,size(imagePixels,1)*size(imagePixels,2)),size(imagePixels,1),size(imagePixels,2));
%     liveWireEngine.setSeed(200,200);
%     returnedPath = liveWireEngine.returnPath(400, 400);
    imshow(mat2gray(imagePixels));
    hold on;
    set(gcf,'position',[10,10,1000,1000]);
    set(gcf,'WindowButtonUpFcn',@mouseLeftClick);  %%LiveWire init and setting points are handled with callbacks
    disp('Callback set');
%     %set(gca,'ydir','reverse');
%     hold on;
%     plot(returnedPath(2,:),returnedPath(1,:),'r');
end