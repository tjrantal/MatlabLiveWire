%%Script mostly written. Perhaps there is an error in indexing. Not all
%%Of the cells get visited. Also, Wherefroms do not form logical paths

%%Check DijsktraheapPorting.java WhereFrom and go from there. Should
%%Have one larger index values + make those match with one larger
%%X and Y indices

function ok = liveWireSetTarget(targetPoint)
    global whereFrom visited parameters;
    parameters.targetPoint = targetPoint;
     
    [returnedPath] = returnPath(parameters.seedPoint,parameters.targetPoint,visited,whereFrom,size(parameters.gradientr,1),size(parameters.gradientr,2));
    parameters.returnedPath = returnedPath;
  

    
    %%%FUNCTIONS%%%%
    
    %Calculate path from seed to x,y
    function [returnedPath] = returnPath(sourcePoint,targetPoint,visited,whereFrom,width,height)
        if visited(targetPoint(1),targetPoint(2)) == 0
           disp('Not there yet');
           ok = 1;
           returnedPath = parameters.returnedPath;
            return; %SetPoint has not finished yet?
           
        end
%         disp('was visited');
        sr = sourcePoint(1);
        sc = sourcePoint(2);
        myr = targetPoint(1);
        myc = targetPoint(2);
        
        length = 0;
%         disp(['C 1 Where from ' num2str(whereFrom(myr,myc)) ' check loop R ' num2str(myr) ' C ' num2str(myc) ' length ' num2str(length)]);
        while myc~=sc || myr~=sr
            length = length+1;
    		nextc = mod(whereFrom(myr,myc),width)+1;
    		nextr = floor(whereFrom(myr,myc)/width)+1;
    		myc = nextc;
    		myr = nextr;
%             disp(['Where from ' num2str(whereFrom(myr,myc)) ' check loop R ' num2str(myr) ' C ' num2str(myc) ' length ' num2str(length)]);
        end
%         disp('check 2');
        mylength = length;
        myr = targetPoint(1);
        myc = targetPoint(2);
        count = 0;
        tempc(1) = myc;
        tempr(1) = myr;
%         disp(['C 2 Where from ' num2str(whereFrom(myr,myc)) ' check loop R ' num2str(myr) ' C ' num2str(myc) ' length ' num2str(length)]);
         while myc~=sc || myr~=sr
            count = count+1;
    		nextc = mod(whereFrom(myr,myc),width)+1;
    		nextr = floor((whereFrom(myr,myc))/width)+1;
            tempr(count+1) = nextr;
            tempc(count+1) = nextc;
    		myc = nextc;
    		myr = nextr;
%             disp(['C 2 Where from ' num2str(whereFrom(myr,myc)) ' check loop R ' num2str(myr) ' C ' num2str(myc) ' length ' num2str(length)]);
         end
        returnedPath= [flipud(tempr);flipud(tempc)];
        ok = 1;
    end
    
 
end