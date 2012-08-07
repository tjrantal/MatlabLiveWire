function testJava
    
    close all;
    clear all;
    javaaddpath('.');
    global lineHandle returnedPath liveWireEngine imagePixels;
    %dicomFileIn = 'IM000~10';
%     dicomFileIn = 'C:\MyTemp\oma\Timon\tyo\SubchondralPilot\karsittu\kh1\18834435';
    dicomFileIn = 'C:\MyTemp\oma\Timon\tyo\SubchondralPilot\livewireData\10022712\18830471';
    imageInfo = dicominfo(dicomFileIn);
    imagePixels = double(dicomread(imageInfo));
    liveWireEngine = javaEngineLiveWire.LiveWireCosts(reshape(imagePixels,1,size(imagePixels,1)*size(imagePixels,2)),size(imagePixels,1),size(imagePixels,2));
    gradientR = liveWireEngine.getGradientR();
    figure;
    subplot(1,2,1);
    imshow(imagePixels,[]);
    subplot(1,2,2);
    imshow(gradientR,[]);
    
    figure
    imshow(mat2gray(imagePixels));
    hold on;
    set(gcf,'position',[10,10,1000,1000]);
    set(gcf,'WindowButtonUpFcn',@mouseLeftClick);  %%LiveWire init and setting points are handled with callbacks
    disp('Callback set');

end