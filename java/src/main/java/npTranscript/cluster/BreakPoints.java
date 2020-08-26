package npTranscript.cluster;

import java.io.IOException;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

class BreakPoints{
 private SparseRealMatrix[][] breakpoints;
 private SparseVector[][] breakSt, breakEnd ;
public void refresh(int seqlen){
	for(int i=0; i<breakSt.length; i++){
		
		for(int j=0; j<breakSt[i].length; i++){
			this.breakpoints[i][j] = new OpenMapRealMatrix(seqlen, seqlen);
			breakSt[i][j].clear();
			breakEnd[i][j].clear();
		}
	}
}
 BreakPoints(int num_sources, int seqlen){
			this.breakpoints = new SparseRealMatrix[num_sources][2];
			this.breakSt = new SparseVector[num_sources][2];
			this.breakEnd = new SparseVector[num_sources][2];
				for(int i=0; i<breakpoints.length; i++){
					this.breakSt[i] = new SparseVector[] {new SparseVector(),new SparseVector()};
					this.breakEnd[i] = new SparseVector[] {new SparseVector(),new SparseVector()};
					breakpoints[i] = new OpenMapRealMatrix[] {new OpenMapRealMatrix(seqlen, seqlen),new OpenMapRealMatrix(seqlen, seqlen)};
				
				}
			
	 }
 void addBreakPoint(int source_index,int i, int prev_position, int position) {
	if(breakpoints[source_index][i]!=null){
		this.breakpoints[source_index][i].addToEntry(prev_position, position, 1);
		this.breakSt[source_index][i].addToEntry(prev_position, 1);
		this.breakEnd[source_index][i].addToEntry(position,  1);
	}
 }
	
		 void printBreakPoints(Outputs o, int chrom_index)  {
			 Runnable run = new Runnable(){
				 
			 public void run(){
			for(int i=0; i<breakpoints.length; i++){
				if(breakpoints[i]!=null && breakpoints[i][0]!=null){
					for(int j=0 ; j<breakpoints[i].length; j++)
					{
				try{
					o.printMatrix(breakpoints[i][j],breakSt[i][j], breakEnd[i][j], chrom_index, i, j);
				}catch(IOException exc){
					exc.printStackTrace();
				}
					
					}
				}
			}
			 }
			 };
			 Outputs.executor.execute(run);
		}
	
	
}