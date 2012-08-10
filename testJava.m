% function digitizedPath = testJava
    
    close all;
    clear all;
    clc;
    javaaddpath('.');
    global lineHandle returnedPath liveWireEngine imagePixels digitizedPath;
    %dicomFileIn = 'IM000~10';
%     dicomFileIn = 'C:\MyTemp\oma\Timon\tyo\SubchondralPilot\karsittu\kh1\18834435';

%     dicomFileIn = 'C:\MyTemp\oma\Timon\tyo\SubchondralPilot\livewireData\10022712\18830471';
    dicomFileIn = 'C:\MyTemp\oma\Timon\tyo\SubchondralPilot\livewireData\10022712\18831575';
    imageInfo = dicominfo(dicomFileIn);
    imagePixels = double(dicomread(imageInfo));
    %Filter with Gaussian
    gaus = fspecial('gaussian',3);
    imagePixels = imfilter(imagePixels,gaus,'replicate');
    liveWireEngine = javaEngineLiveWire.LiveWireCosts(imagePixels);
%     gradientr = liveWireEngine.gradientr;
%     figure;
%     subplot(1,3,1);
%     imshow(imagePixels,[]);
%     subplot(1,3,2);
%     imshow(gradientr,[]);
%     laplacian = liveWireEngine.laplacian();
%     subplot(1,3,3);
%     imshow(laplacian,[]);
%     keyboard
    figure
    imshow(mat2gray(imagePixels));
    hold on;
    set(gcf,'position',[10,10,1000,1000]);
    set(gcf,'WindowButtonUpFcn',@mouseLeftClick);  %%LiveWire init and setting points are handled with callbacks
    disp('Callback set');

% end