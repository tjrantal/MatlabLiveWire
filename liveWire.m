
function parameters = liveWire(dicomFileIn);
    global lineHandle parameters imagePixels whereFrom seedRunning keepGoing;
    imageInfo = dicominfo(dicomFileIn);
    imagePixels = double(dicomread(imageInfo));
    [parameters] = initGradient(imagePixels);
    parameters.gw = 0.43;
    parameters.dw = 0.13;
    parameters.ew = 0.0;
    parameters.pw = 30;
    figure;
    %pcolor(parameters.gradientr)
    %imshow(parameters.gradientr)
    disp('Init ready');
    imshow(mat2gray(imagePixels));
    set(gcf,'position',[10,10,1000,1000]);
    %set(gca,'ydir','reverse');
    hold on;
    disp('Callback set');
    set(gcf,'WindowButtonUpFcn',@mouseLeftClick);  %%LiveWire init and setting points are handled with callbacks
end
     