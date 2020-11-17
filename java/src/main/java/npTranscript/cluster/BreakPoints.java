package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;

import japsa.seq.Sequence;
import npTranscript.NW.AlignmentParameters;
import npTranscript.NW.AlignmentResult;
import npTranscript.NW.NeedlemanWunsch;

class BreakPoints{
	public static boolean writeScores = true;
	
 private OpenMapRealMatrix[][] breakpoints;
 final Sequence seq;
 List<Sequence> genomes;
 final int num_sources;
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

SortedSet<Integer> getPositions(SparseVector[][] starts,  int lev){
	TreeSet<Integer> s = new TreeSet<Integer>();
	for(int i=0; i<num_sources; i++){
		s.addAll(starts[i][lev].m.keySet());
	}
	return s;
}
public boolean contains(int start, int end, int lev){
	for(int i=0; i<num_sources; i++){
		if(this.breakpoints[i][lev].getEntry(start, end)>0){
			return true;
		}
	}
	return false;
}
public SparseRealMatrix calculateScores(int k, int[] left, int[] right, Sequence left_sequence, Sequence right_sequence){
	SortedSet<Integer> starts = getPositions(this.breakSt,k);
	SortedSet<Integer> ends = getPositions(this.breakEnd,k);
	AlignmentParameters alignP=new AlignmentParameters(2,1);
	 int seqlen = seq.length();
	SparseRealMatrix scores  = new OpenMapRealMatrix(seqlen, seqlen);
	for(Iterator<Integer> it = starts.iterator(); it.hasNext();){
		Integer st = it.next();
		Sequence leftSeq = left_sequence.subSequence(st + left[0], st+left[1]); // not sure if we should substract one
		for(Iterator<Integer> it1= ends.iterator(); it1.hasNext();){
			Integer en = it1.next();
			//if(contains(st,en,k)){
				Sequence rightSeq = right_sequence.subSequence(en + right[0], en+right[1]);
				AlignmentResult aligns = NeedlemanWunsch.computeNWAlignment(leftSeq.toString(), rightSeq.toString(),   alignP);
				double perc = 100.0/(double) Math.exp(aligns.getTotalCost());//aligns.getMatchPerc()*100;
				scores.setEntry(st, en, perc);;
		//	}
			//aligns.getTotalCost();
		}
	}
	return scores;
}
//this is in a 1 based system, so need add an extra position
 BreakPoints(int num_sources, Sequence seq, ArrayList<Sequence> genomes){
	 this.seq = seq;
	 this.genomes = genomes;
	 this.num_sources = num_sources;
	 int seqlen = seq.length()+1;
	 
	 int num = genomes.size()+1;
	 for(int i=0; i<genomes.size(); i++){
		 chr_names.put(genomes.get(i).getName(), i);
	 }
			this.breakpoints = new OpenMapRealMatrix[num_sources][num];
			this.breakSt = new SparseVector[num_sources][num];
			this.breakEnd = new SparseVector[num_sources][num];
				for(int i=0; i<breakpoints.length; i++){
					this.breakSt[i] = new SparseVector[num] ;
					this.breakEnd[i] = new SparseVector[num];
					breakpoints[i] = new OpenMapRealMatrix[num];
					for(int j=0; j<num; j++){
						this.breakSt[i][j] = new SparseVector();
						this.breakEnd[i][j] = new SparseVector();
						this.breakpoints[i][j] = new OpenMapRealMatrix(seqlen, seqlen);
					}
				}
			
	 }
 void addBreakPoint(int source_index,int i, int prev_position, int position) {
	if(breakpoints[source_index][i]!=null){
		this.breakpoints[source_index][i].addToEntry(prev_position, position, 1);
		this.breakSt[source_index][i].addToEntry(prev_position, 1);
		this.breakEnd[source_index][i].addToEntry(position,  1);
	}
 }
 final Map<String, Integer> chr_names = new HashMap<String, Integer>();
 /** this is for adding fusion breakpoints */
 void addBreakPoint(int source_index,String str, int prev_position, int position) {
	 int i = 1+chr_names.get(str);
		if(breakpoints[source_index][i]!=null){
			this.breakpoints[source_index][i].addToEntry(prev_position, position, 1);
			this.breakSt[source_index][i].addToEntry(prev_position, 1);
			this.breakEnd[source_index][i].addToEntry(position,  1);
		}
	 }
 

	
		 void printBreakPoints(Outputs o, int chrom_index)  {
			 Runnable run = new Runnable(){
			 public void run(){
				 int num = breakpoints[0].length;
				 if(breakpoints[0]!=null){
					 for(int j=0 ; j<num; j++){
						 Sequence left_seq = j<2 ? seq : genomes.get(j-1);
						 Sequence right_seq = seq;
					 SparseRealMatrix scores55 =  calculateScores(j, new int [] {-10,0}, new int[] {-10,0}, left_seq, right_seq);
					 SparseRealMatrix scores33 =  calculateScores(j, new int [] {0,10}, new int[] {0,10}, left_seq, right_seq);
					 SparseRealMatrix scores53 =  calculateScores(j, new int [] {-10,0}, new int[] {0,10}, left_seq, right_seq);
					 SparseRealMatrix scores35 =  calculateScores(j, new int [] {0,10}, new int[] {-10,0}, left_seq, right_seq);
					 SparseRealMatrix scoresMid =  calculateScores(j, new int [] {-5,5}, new int[] {-5,5}, left_seq, right_seq);

					 
					 for(int i=0; i<breakpoints.length; i++){
						if(breakpoints[i]!=null && breakpoints[i][j]!=null){
						try{
							o.printMatrix(breakpoints[i][j],breakSt[i][j], breakEnd[i][j],  chrom_index, i, j, "");
							o.printMatrix(scores55,breakSt[i][j], breakEnd[i][j],  chrom_index, i, j, "scores55/");
							o.printMatrix(scores33,breakSt[i][j], breakEnd[i][j],  chrom_index, i, j, "scores33/");
							o.printMatrix(scores53,breakSt[i][j], breakEnd[i][j],  chrom_index, i, j, "scores53/");
							o.printMatrix(scores35,breakSt[i][j], breakEnd[i][j],  chrom_index, i, j, "scores35/");
							o.printMatrix(scores35,breakSt[i][j], breakEnd[i][j],  chrom_index, i, j, "scoresMid/");


						}catch(IOException exc){
							exc.printStackTrace();
						}
						
						}
					}
				 }
			}
					o.breakPW.close();
			 }
			 };
			 Outputs.h5writer.execute(run);
		}
	
	
}