    %Init gradients, used for cost function
    function [parameters] =initGradient(imagePixels)
        %Sobel filter kernels
        Gx = [-1,0,1;
              -2,0,2;
              -1, 0, 1;];
          %Y increases from top to bottom -> convolution kernel is flipped
          Gy = [1,2,1;
              0,0,0;
              -1, -2, -1;];
          parameters.gradientx = conv2(imagePixels,Gx,'same');
          parameters.gradienty = conv2(imagePixels,Gy,'same');
          
          parameters.gradientr = sqrt(parameters.gradientx.^2+parameters.gradienty.^2);
          parameters.grmin = min(min(parameters.gradientr));
          parameters.grmax = max(max(parameters.gradientr));         
    end