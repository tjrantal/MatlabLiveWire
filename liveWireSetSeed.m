%%Script mostly written. Perhaps there is an error in indexing. Not all
%%Of the cells get visited. Also, Wherefroms do not form logical paths

%%Check DijsktraheapPorting.java WhereFrom and go from there. Should
%%Have one larger index values + make those match with one larger
%%X and Y indices

function [parametersReturn] = liveWireSetSeed
    %these are default parameter values

    global whereFrom visited keepGoing seedRunning parameters;    
    pixelCosts   = struct([]); 
    seedPoint   =parameters.seedPoint;
    setPoint(parameters.gradientx,parameters.gradienty,parameters.gradientr,parameters.grmin,parameters.grmax,seedPoint,pixelCosts,parameters.gw,parameters.dw,parameters.ew,parameters.pw);
    parameters.visited      = visited;
    parameters.whereFrom    = whereFrom;
    
    %%%FUNCTIONS%%%%
    

    %Set seed point and calculate costs to all points
    function setPoint(gradientx,gradienty,gradientr,grmin,grmax,seedPoint,pixelCosts,gw,dw,ew,pw)
        r = seedPoint(1);
        c = seedPoint(2);
        sr = seedPoint(1);
        sc = seedPoint(2);
        width = size(gradientx,1);
        height = size(gradientx,2);
        visited = zeros(width,height,'uint8');   %For maintaining pixels gone through in Dijsktra
        whereFrom = zeros(width,height); %For source in Dijkstra
        
        visited(r,c) = 1;
        whereFrom(r,c) = width*(r-1)+c-1;
        disp('Into setting seed');
        [visited,pixelCosts] = updateCosts(width,height,r,c,0,visited,pixelCosts,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
%          tempFig = figure('position',[10 10 1000 1000]);
%          imageHandle = imshow(mat2gray(visited));
        while(length(pixelCosts) > 0 && keepGoing == 1)
            %pixelCosts
            nextPixel = peekQueue(pixelCosts);
            nextIndex = nextPixel.index;
            nextc = mod(nextIndex,width)+1;
			nextr = floor((nextIndex)/width)+1;
%             disp([ 'Visiting index ' num2str(nextIndex) ' X ' num2str(nextc) ' Y ' num2str(nextr) ' from '  num2str(nextPixel.whereFrom)]);
            whereFrom(nextr,nextc) =nextPixel.whereFrom;
            [visited,pixelCosts] = updateCosts(width,height,nextr, nextc, nextPixel.cost,visited,pixelCosts,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
            %removes pixels that are already visited and went to the queue
            goOn = 1;
%              set(imageHandle,'CData',mat2gray(visited));
%              drawnow
            while goOn ==1 && keepGoing == 1
                nextPixel = peekQueue(pixelCosts);
                if( ~isstruct(nextPixel))
                    goOn = 0;
                else
                    nextIndex = nextPixel.index;
                    nextc = mod(nextIndex,width)+1;
                    nextr = floor((nextIndex)/width)+1;
                    %disp(['Testing X ' num2str(nextX) ' Y ' num2str(nextY) ' from '  num2str(nextPixel.whereFrom)]);
                    if(visited(nextr,nextc)==0)
                        goOn = 0;
                    else
                        pixelCosts = pollQueue(pixelCosts);  
                        %disp(['Removing X ' num2str(nextX) ' Y ' num2str(nextY) ' from '  num2str(nextPixel.whereFrom)]);
                    end
                end

            end

        end
        disp('Stopped calculating costs');
        seedRunning = 0;
%         close(tempFig);
        while length(pixelCosts) > 0
            pixelCosts = pollQueue(pixelCosts); 
        end       
    end
    
    %Remove smallest value from queue
    function [pixelCosts] = pollQueue(pixelCosts)
        if length(pixelCosts) > 0
           [val ind] = sort([pixelCosts.cost]);
           pixelCosts(ind(1)) = [];
        end
    end

    function [nextPixel] = peekQueue(pixelCosts)
        if length(pixelCosts) > 0
           [val ind] = sort([pixelCosts.cost]);
           nextPixel = pixelCosts(ind(1));
        else
            nextPixel = 0;
        end
    end
    %Calculate costs
    %Added constraint that only pixels less than 10 pixels away are
    %considered
    function [visited,pixelCosts] = updateCosts(width,height,r,c,mycost,visited,pixelCosts,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw)
        maxDist = 20;
        pixelCosts = pollQueue(pixelCosts);
        visited(r,c) = 1;
        %upper right (r-1,c+1)
        if (c < width ) && (r > 1) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (r-1-1)*width+c+1-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r-1,c+1,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end
        
        %upper left (r-1,c-1)
        if (c > 1) && (r > 1) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c-1)+(r-1-1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r-1,c-1,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end

		%down right (r+1,c+1)
        if (c < width) && (r < height) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c+1)+(r-1+1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r+1,c+1,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end
        
		%down left (r+1,c-1)
        if (c > 1) && (r < height) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c-1)+(r-1+1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r+1,c-1,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end
		
		%update left cost (r,c-1)
        if (c > 1) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c-1)+(r-1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r,c-1,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end

        %update right cost (r,c+1)
        if (c < width) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c+1)+(r-1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r,c+1,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end
		
		%update up cost (r-1,c)
        if (r > 1) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c)+(r-1-1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r-1,c,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end
        
		%update down cost (r+1,c)
        if (r < height) && (sqrt((parameters.seedPoint(1)-r)^2+(parameters.seedPoint(2)-c)^2)<=maxDist)
           indeksi = length(pixelCosts)+1;
           pixelCosts(indeksi).index = (c)+(r+1-1)*width-1;
           pixelCosts(indeksi).cost = mycost+edgeCost(r,c,r+1,c,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw);
           pixelCosts(indeksi).whereFrom = (r-1)*width+c-1;
        end
    end

    %Edge cost
    function returnCost = edgeCost(sr,sc,dr,dc,gradientx,gradienty,gradientr,grmin,grmax,gw,dw,ew,pw)
        %fg is the Gradient Magnitude
		%we are dividing by sqrt(2) so that the value won't pass 1
		%as is stated in United Snakes formule 36
        fg = (1.0/sqrt(2)*sqrt( (dc-sc)*(dc-sc) + (dr-sr)*(dr-sr))* (1 - ((gradientr(dr,dc)-grmin)/(grmax-grmin))));

		if grmin==grmax
			fg= (1.0/sqrt(2)*sqrt( (dc-sc)*(dc-sc) + (dr-sr)*(dr-sr)));
        end
        
        temp = (gradientr(dr,dc)-grmin)/(grmax-grmin);
		fe = exp(-pw*temp)*fg;
        
        %Gradient direction
        GradVector = [gradientx(sr,sc), gradienty(sr,sc)];  %Gradient direction vector
        Dp = GradVector/norm(GradVector);                   %Gradient direction Unit vector
        DpN = [Dp(2) -Dp(1)];                               %Normal vector to Dp
        
        %Lpq is the normalized biderectional link between pixels p and q (United Snakes formule 38)
		p = [sr,sc];
		q = [dr,dc];
		
		%remember that y is upside down... we need to invert the y term in pq
		
		myPQ = [dr - sr, sc - dc];
		
		
		if(dot(DpN,myPQ ) >=0 )
			Lpq = myPQ/norm(myPQ); %(q-p)/||p-q|| 
        else	    
			Lpq = [sc-dc,dr-sr]/norm([sc-dc,dr-sr]); % (p-q)/||p-q||
        end
		
		dppq = dot(DpN,Lpq);

		GradVectorq = [gradientx(dr,dc),gradienty(dr,dc)];
		Dq = GradVectorq/norm(GradVectorq);
		DqN = [Dq(2) -Dq(1)];
		dqpq = dot(Lpq,DqN);

		%United Snakes formule 37
		%I have found a problem here... 
		%When the gradient is near zero in the place, when getting the unit vector, 
		%it becomes NaN, because we are dividing something for about zero... 
		%but the value of fd should be as high as possible (we are over an edge)
		%so, when asking for unit vectors, we check for components x and y. If they are
		%too small, we return a vector (1,0)
		if (abs(GradVector(1))<0.0000001)&&(abs(GradVector(2))<0.0000001)
			dppq = 0;
        end
		if(abs(GradVectorq(1))<0.0000001)&&(abs(GradVectorq(2))<0.0000001)
			dqpq = 0;
        end
		
		fd = 2.0/(3*pi)*(acos(dppq)+acos(dqpq));	    
		%return Math.abs(imagePixels[toIndex(dx,dy)]-imagePixels[toIndex(sc,sr)]);
		returnCost = ew*fe+gw*fg+dw*fd;%+0.2*Math.sqrt( (dx-sc)*(dx-sc) + (dy-sr)*(dy-sr));
        
    end
end