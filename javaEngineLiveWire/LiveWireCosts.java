/*
	Modified by Timo Rantalainen from IvusSnakes ImageJ plugin 
	A Class to calculate LiveWire paths for Matlab. Matlab didn'
	provide reasonable speed -> replaced Matlab with java to
	speed up things...
	
	Changed the implementation back to the one suggested in Barret & Mortensen 1997.
	Interactive live-wire boundary extraction. Medical Image Analysis (1996/7) volume 1, number 4, pp 331–341.
	
	Use:
		Create object with pixel data and image width and size:
			object = LiveWireCosts(pixelMatrix,width,height);
		set seed point:
			object.setSeed(x,y);
*/

package javaEngineLiveWire;

import java.util.PriorityQueue;


public class LiveWireCosts implements Runnable{
	
	
	
    double[] imagePixels; //stores Pixels from original image
    int[] imageCosts; //stores Costs for every pixel
    PriorityQueue<PixelNode> pixelCosts;
    double[] gradientx; //stores image gradient modulus 
    double[] gradienty; //stores image gradient modulus 
    //it is oriented: X = LEFT TO RIGHT
    //                Y = UP   TO DOWN
    double[] gradientr; //stores image gradient RESULTANT modulus 
	double[] laplacian;
    double grmin;//gradient global minimum
    double grmax;//gradient global maximum

    int[] whereFrom;  //stores where from path started
    boolean[] visited; //stores whether the node was marked or not
    int width;
    int height;
    int sx,sy; //seed x and seed y, weight zero for this point

    private int tx,ty;//thread x and y passed as parameters

    Thread myThread;
    boolean myThreadRuns;//flag for thread state
    
    private double gw;//Gradient Magnitude Weight - set by setGWeight
    private double dw;//Gradient Direction Weight - set by setDWeight
	private double zw;//Exponential        Weight - set by setEWeight
    private double ew;//Exponential        Weight - set by setEWeight
    private double pw;//Exponential Potence Weight - set by setPWeight
    
    static int INF = 0x7FFFFFFF; //maximum integer

    //converts x, y coordinates to vector index
    private int toIndex(int x,int y){
		return (y*width+x);
    }

    //initializes gradient vector
    private void initGradient(){
		gradientx = new double[height*width];
		gradienty = new double[height*width];
		gradientr = new double[height*width];
		//Using sobel
		//for gx convolutes the following matrix
		//   
		//     |-1 0 1|
		//Gx = |-2 0 2|
		//     |-1 0 1|


		for(int i=0;i<width;i++){
			for(int j=0;j<height;j++){
			if((i>0)&&(i<width-1)&&(j>0)&&(j<height-1)){
				gradientx[toIndex(i,j)] =
				-1*(imagePixels[toIndex(i-1,j-1)]) +1*(imagePixels[toIndex(i+1,j-1)])
				-2*(imagePixels[toIndex(i-1,j  )]) +2*(imagePixels[toIndex(i+1,j  )])
				-1*(imagePixels[toIndex(i-1,j+1)]) +1*(imagePixels[toIndex(i+1,j+1)]);
			}
			}
		}

		//for gy convolutes the following matrix (remember y is zero at the top!)
		//
		//     |+1 +2 +1| 
		//Gy = | 0  0  0|
		//     |-1 -2 -1|
		//
		for(int i=0;i<width;i++){
			for(int j=0;j<height;j++){
			if((i>0)&&(i<width-1)&&(j>0)&&(j<height-1)){
				gradienty[toIndex(i,j)] = 
				+1*(imagePixels[toIndex(i-1,j-1)]) -1*(imagePixels[toIndex(i-1,j+1)])
				+2*(imagePixels[toIndex(i  ,j-1)]) -2*(imagePixels[toIndex(i  ,j+1)])
				+1*(imagePixels[toIndex(i+1,j-1)]) -1*(imagePixels[toIndex(i+1,j+1)]);
			}
			}
		}
		for(int i=0;i<width;i++){
			for(int j=0;j<height;j++){
			if((i>0)&&(i<width-1)&&(j>0)&&(j<height-1)){
				//Math.hypot returns sqrt(x^2 +y^2) without intermediate overflow or underflow.
				
				gradientr[toIndex(i,j)] = Math.sqrt( gradientx[toIndex(i,j)]*gradientx[toIndex(i,j)]+
								 gradienty[toIndex(i,j)]*gradienty[toIndex(i,j)]);
				
			}
			}
		}

		grmin = gradientr[0];
		grmax = gradientr[0];
		for(int i=0;i< height*width;i++){
			if(gradientr[i]<grmin) grmin=gradientr[i];
			if(gradientr[i]>grmax) grmax=gradientr[i];
		}
    }

    public double[][] getGradientR(){
		double[][] gradientR = new double[width][height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				gradientR[j][i] = gradientr[i*width+j];
			}
		}
		return gradientR;
    }
	
	public double[][] getLaplacian(){
		double[][] laplacianR = new double[width][height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				laplacianR[j][i] = laplacian[i*width+j];
			}
		}
		return laplacianR;
    }

	    //initializes gradient vector
    private void initLaplacian(){
		laplacian = new double[height*width];
		
		//Using finite differences
		//convolute with
		//   
		//     |0	1 0|
		//Gx = |1	-4 1|
		//     |0	1 0|
		double[][] laplacianKernel = {{0,1,0},
									{1,-4,1},
									{0,1,0}};

		for(int i=1;i<width-1;i++){
			for(int j=1;j<height-1;j++){
				for (int j2=-1;j2<=1;++j2){
					for (int i2=-1;i2<=1;++i2){
						laplacian[toIndex(i,j)] += imagePixels[toIndex(i+i2,j+j2)]*laplacianKernel[i2+1][j2+1];
					}
				}
			}
		}

		/*Search for zero crossing to binarize the result*/
		double[] tempLap = new double[laplacian.length];
		int[][] neighbourhood = new int[8][2];	//8 connected neighbourhood
		for(int i=1;i<width-1;i++){
			for(int j=1;j<height-1;j++){
				tempLap[toIndex(i,j)] = 1;	
				if (laplacian[toIndex(i,j)] == 0){	/*No need to check neighbours*/
					tempLap[toIndex(i,j)] = 0;	
				}else{	/*Check neighbours*/
					//Check 8-connected neighbour
					for (int th = 0; th<8;++th){
						neighbourhood[th][0] = i+(int) Math.round(Math.cos((double) th));
						neighbourhood[th][1] = j+(int) Math.round(Math.sin((double) th));
					}
					int[] centre = {i,j};
					tempLap = checkNeighbours(tempLap,neighbourhood,centre);
					
				}
			}
		}
		
		/*OverWrite Laplacian*/
		for(int i=0;i<width;i++){
			for(int j=0;j<height;j++){
				laplacian[toIndex(i,j)]  = tempLap[toIndex(i,j)];	
			}
		}
    }
	
		/*Update pixel queue*/
	protected double[] checkNeighbours(double[] tempLap, int[][] neighbourhood,int[] centre){
		int[] coordinates;
		double[] lbpHist;
        for (int r = 0;r<neighbourhood.length;++r){
			coordinates = neighbourhood[r];
            if (Math.signum(laplacian[toIndex(coordinates[0],coordinates[1])]) != Math.signum(laplacian[toIndex(centre[0],centre[1])])){ /*Signs differ, mark border*/
				if (Math.abs(laplacian[toIndex(centre[0],centre[1])]) < Math.abs(laplacian[toIndex(coordinates[0],coordinates[1])])){
					tempLap[toIndex(centre[0],centre[1])] = 0;
			    }else{
					tempLap[toIndex(coordinates[0],coordinates[1])] = 0;
				}
            }
        }
		return tempLap;
	}

	
    //initializes Dijkstra with the image
    public LiveWireCosts(double[] image,int x, int y){
 
		//initializes weights for edge cost taken from Barret 1997
    	//these are default values
    	gw = 0.43;
		zw = 0.43;
    	dw = 0.13;
    	
		//initializes all other matrices
		imagePixels = new double[x*y];
		pixelCosts = new PriorityQueue<PixelNode>();
		whereFrom   = new int [x*y];
		visited     = new boolean[x*y];
		width  = x;
		height = y;
				
		//copy image matrice
		for(int j=0;j<y;j++){
			for(int i=0;i<x;i++){	    	
			imagePixels[j*x+i] = 
				image[j*x+i];		
			//imageCosts [j*x+i] = INF;
			visited    [j*x+i] = false;
			}
		}
		initGradient();
		initLaplacian();
		//inits the thread
		myThread = new Thread(this);

    }    
    //returns de cost of going from sx,sy to dx,dy
    private double edgeCost(int sx,int sy,int dx,int dy){
		//fg is the Gradient Magnitude

		//we are dividing by sqrt(2) so that the value won't pass 1
		//as is stated in United Snakes formule 36
		
		double fg = (1.0/Math.sqrt(2)*Math.sqrt( (dx-sx)*(dx-sx) + (dy-sy)*(dy-sy))* 
			(1 - ((gradientr[toIndex(dx,dy)]-grmin)/(grmax-grmin))));

		if(grmin==grmax) 
			fg= (1.0/Math.sqrt(2)*Math.sqrt( (dx-sx)*(dx-sx) + (dy-sy)*(dy-sy)));


		//this parameter is an attempt to find edges in IVUS images	
		double x = (gradientr[toIndex(dx,dy)]-grmin)/(grmax-grmin);
		double fe = Math.exp(-pw*x)*fg;


		//	System.out.println("Fg " + fg +" gradientr " + gradientr[toIndex(dx,dy)] + " grmin " + grmin + " grmax " + grmax);
		//fd id the Gradient Direction

		//CHECK ME OUT: The part of fd has not been much tested
		//if someone wishes to spend some time testing it, it'd be a great idea

		//Dp is the unit vector of the gradient direction at pixel p (sx,sy)
		//this is defined near formule 37 in United Snakes
		Vector2d GradVector = new Vector2d(gradientx[toIndex(sx,sy)],gradienty[toIndex(sx,sy)]);
		Vector2d Dp = GradVector.getUnit();
		//DpN is the normal vector to Dp
		Vector2d DpN = new Vector2d(Dp.getNormal());
		//Lpq is the normalized biderectional link between pixels p and q (United Snakes formule 38)
		Vector2d Lpq = new Vector2d();
		Vector2d p = new Vector2d(sx,sy);
		Vector2d q = new Vector2d(dx,dy);
		
		//remember that y is upside down... we need to invert the y term in pq
		
		Vector2d myPQ = new Vector2d ( dx - sx, sy - dy);
		
		
		if(DpN.dotProduct( myPQ ) >=0 ){
			Lpq = myPQ.getUnit(); // (q-p)/||p-q|| 
		}
		else{	    
			Lpq = myPQ.getUnit(); // (p-q)/||p-q||
		}
		//dppq = DpN . Lpq 
		double dppq = DpN.dotProduct(Lpq);
		//dqpq = Lpq . DpN
		Vector2d GradVectorq = new Vector2d(gradientx[toIndex(dx,dy)],gradienty[toIndex(dx,dy)]);
		Vector2d Dq = GradVectorq.getUnit();
		Vector2d DqN = new Vector2d(Dq.getNormal());
		double dqpq = Lpq.dotProduct(DqN);

		//United Snakes formule 37
		//I have found a problem here... 
		//When the gradient is near zero in the place, when getting the unit vector, 
		//it becomes NaN, because we are dividing something for about zero... 
		//but the value of fd should be as high as possible (we are over an edge)
		//so, when asking for unit vectors, we check for components x and y. If they are
		//too small, we return a vector (1,0)
		if((Math.abs(GradVector.getX())<0.0000001)&&(Math.abs(GradVector.getY())<0.0000001))
			dppq = 0;
		if((Math.abs(GradVectorq.getX())<0.0000001)&&(Math.abs(GradVectorq.getY())<0.0000001))
			dqpq = 0;
			
		
		double fd = 2.0/(3*Math.PI)*(Math.acos(dppq)+Math.acos(dqpq));	    

		/*	System.out.println("Fd "+ fd + " Fg " + fg + " acos dppq " + Math.acos(dppq) + " acos dqpq " + Math.acos(dqpq)
				   + " Dp (" + Dp.getX() + "," + Dp.getY() +") DpN (" + DpN.getX() +"," + DpN.getY() + ")"  
				   + " p (" +p.getX()+ ","+p.getY() + ") q(" +q.getX()+","+q.getY()+") " 
				   + "Gradp(" + gradientx[toIndex(sx,sy)]+ ","+gradienty[toIndex(sx,sy)]+ ") " 
				   + "Gradq(" + gradientx[toIndex(dx,dy)]+ ","+gradienty[toIndex(dx,dy)]+ ")");*/

		//return Math.abs(imagePixels[toIndex(dx,dy)]-imagePixels[toIndex(sx,sy)]);
		return ew*fe+gw*fg+dw*fd;//+0.2*Math.sqrt( (dx-sx)*(dx-sx) + (dy-sy)*(dy-sy));
		
    }
    //sets weights for fg and fd
    public void setGWeight(double pgw){
    	gw = pgw;
    }
    public void setDWeight(double pdw){
    	dw = pdw;
    }
    public void setEWeight(double pew){
    	ew = pew;
    }
    public void setPWeight(double ppw){
    	pw = ppw;
    }
    
    //updates Costs and Paths for a given point
    //actuates over 8 directions N, NE, E, SE, S, SW, W, NW
    private void updateCosts(int x,int y,double mycost){

	visited[toIndex(x,y)] = true;
	pixelCosts.poll();

//	mycost = mycost;
	//upper right
	if((x< width-1)&&(y>0)){
	    pixelCosts.add(new PixelNode(toIndex(x+1,y-1), mycost+edgeCost(x,y,x+1,y-1),toIndex(x,y)));	    
	}
	//upper left
	if((x>0)&&(y>0)){
	    pixelCosts.add(new PixelNode(toIndex(x-1,y-1), mycost+edgeCost(x,y,x-1,y-1),toIndex(x,y)));	    
	}
	//down right
	if((x< width-1)&&(y<height-1)){
	    pixelCosts.add(new PixelNode(toIndex(x+1,y+1), mycost+edgeCost(x,y,x+1,y+1),toIndex(x,y)));	    
	}
	//down left
	if((x>0)&&(y<height-1)){
	    pixelCosts.add(new PixelNode(toIndex(x-1,y+1), mycost+edgeCost(x,y,x-1,y+1),toIndex(x,y)));	    
	}

	//update left cost
	if(x>0){
	    pixelCosts.add(new PixelNode(toIndex(x-1,y), mycost+edgeCost(x,y,x-1,y),toIndex(x,y)));	    
	}
	//update right cost
	if(x<width-1){
	    pixelCosts.add(new PixelNode(toIndex(x+1,y), mycost+edgeCost(x,y,x+1,y),toIndex(x,y)));
	}
	
	//update up cost
	if(y>0){
	    pixelCosts.add(new PixelNode(toIndex(x,y-1), mycost+edgeCost(x,y,x,y-1),toIndex(x,y)));
	}
	    //update down cost
	if(y<height-1){
	    pixelCosts.add(new PixelNode(toIndex(x,y+1), mycost+edgeCost(x,y,x,y+1),toIndex(x,y)));
	}
    }

    // returns index pointing to next node to be visited
    // It is defined as the minimum cost not yet visited
    // returns -1 if no node is available
    private int findNext(){
	int min = INF;
	int ans = -1;
	for(int y=0;y<height;y++){
	    for(int x=0;x<width;x++){
		if( ( visited   [toIndex(x,y)] == false) &&
		    ( imageCosts[toIndex(x,y)] <  min) ){
		    min = imageCosts[toIndex(x,y)];
		    ans = toIndex(x,y);
		}
	    }
	}
	return ans;
    }

    public int[][] returnPath(int x, int y){
	//returns the path given mouse position

    	int[] tempx = new int[width*height];
		int[] tempy = new int[width*height];
    	
    	if(visited[toIndex(x,y)]==false){
    		//attempt to get path before creating it 
    		//this might occur because of the thread
    		return null;
    	}
    	int length =0;
    	int myx = x;
    	int myy = y;
    	int nextx;
    	int nexty;
    	do{ //while we haven't found the seed	
    		length++;
    		nextx = whereFrom[toIndex(myx,myy)]%width;
    		nexty = whereFrom[toIndex(myx,myy)]/width;
    		myx = nextx;
    		myy = nexty;
    		
    	}while (!((myx==sx)&&(myy==sy)));
    	
    	
    	
    	//add points to vector
    	myx=x;
    	myy=y;
    	int count=0;
    	tempx[0]=myx;//add last points
    	tempy[0]=myy; 
	//	System.out.println("Caminho ");
    	do{ //while we haven't found the seed	    	
    		nextx = whereFrom[toIndex(myx,myy)]%width;
    		nexty = whereFrom[toIndex(myx,myy)]/width;
		//System.out.println("("+nextx+","+nexty+")");
    		
    		count++;
    		tempx[count]=nextx;
    		tempy[count]=nexty; 
    		
    		myx = nextx;
    		myy = nexty;
    		
    	}while (!((myx==sx)&&(myy==sy)));
    	//path is from last point to first
    	//we need to invert it
	//	System.out.println("Caminho ");
		int[][] pathToReturn = new int[2][count+1];
    	for(int i=0;i<=count;i++){

    		pathToReturn[0][i]= tempx[count-i];
    		pathToReturn[1][i]=tempy[count-i];    	
		//System.out.println("( "+vx[i] + " , " + vy[i] + " )");
    	}		    	
    	
    	return pathToReturn;    		    	
    	
    }

    //set point to start Dijkstra
    public void setSeed(int x, int y){
	    	
    	myThreadRuns=false;    	    	
    	try {
			myThread.join();
		} catch (InterruptedException e) {
			System.out.println("Bogus Exception");
			e.printStackTrace();			
		}
    	
    	tx = x;	
    	ty = y;	
    	myThreadRuns = true;
    	myThread = new Thread(this);
		myThread.start();		

    }    
    public void run(){
		//runs set point in parallel
		int x = tx;
		int y = ty;

		int nextIndex;
		int nextX;
		int nextY;
		sx = x;
		sy = y;
		
		for(int i=0;i<height*width;i++){
			//		imageCosts[i]  = INF;
			visited[i]  = false;
		}
		
		
		visited   [toIndex(x,y)] = true; //mark as visited
		//	imageCosts[toIndex(x,y)] = 0; //sets initial point with zero cost
		whereFrom [toIndex(x,y)] = toIndex(x,y);

		

		//update costs
		updateCosts(x,y,0);
		//	nextIndex = findNext();
		//	nextX = nextIndex%width;
		//	nextY = nextIndex/width;
		int debugcount = 0;

		
		while((pixelCosts.peek()!=null)&&(myThreadRuns)){
			//	    System.out.println("Debug count " + debugcount++);
			
			
			
			nextIndex = ((PixelNode)pixelCosts.peek()).getIndex();
			nextX = nextIndex%width;
			nextY = nextIndex/width;
			
			whereFrom[nextIndex] =((PixelNode)pixelCosts.peek()).getWhereFrom();
			
			updateCosts(nextX, nextY, ((PixelNode) pixelCosts.peek()).getDistance());

			//removes pixels that are already visited and went to the queue
			while(true){
			if( pixelCosts.peek() == null )
				break;
			if(visited[ ((PixelNode)pixelCosts.peek()).getIndex() ]==false)
				break;
			pixelCosts.poll();
			}
		}
		while(pixelCosts.peek()!=null)
			pixelCosts.poll();
		
		//	System.out.println("Point set.......");
		/*
		System.out.println("");
		for(int j=0;j<height;j++){
			for(int i=0;i<width;i++){
			System.out.print(imageCosts[j*width+i]+ " ");
			}
			System.out.println("");
			}
		
		System.out.println("Caminhos");
		for(int j=0;j<height;j++){
			for(int i=0;i<width;i++){
			System.out.print("( " + whereFrom[j*width+i]%width + ", " + whereFrom[j*width+i]/height + ") ");
			}
			System.out.println("");
			}*/
		
    }

	public int getTx() {
		return tx;
	}

	public int getTy() {
		return ty;
	}


		    	


}


//this class was created to store pixel nodes in a Priority Queue
//so that Dijkstra can run on O(n log n)
//The interface Comparable is required so that the Java class PriorityQueue
//could be used
//we could not use a standard PriorityQueue<Integer> cause it would 
//

class PixelNode implements Comparable<PixelNode> {
    private int myIndex;
    private double myDistance;
    private int whereFrom;
    public PixelNode(int index, double distance, int whereFrom){
	myIndex = index;
	myDistance = distance;
	this.whereFrom = whereFrom;
    }
    public double getDistance(){
	return myDistance;
    }
    public int getIndex(){
	return myIndex;
    }
    public int getWhereFrom(){
	return whereFrom;
    }

    public int compareTo(PixelNode other){
	if( myDistance < other.getDistance()) 
	    return -1;
	else if( myDistance > other.getDistance()) 
	    return +1;
	else 
	    return 0;
	//	return (int)((myDistance - other.getDistance()));//plus 0.5 to round
    }

}


