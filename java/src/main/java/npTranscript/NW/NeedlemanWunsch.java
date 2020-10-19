package npTranscript.NW;

import java.io.PrintWriter;

public class NeedlemanWunsch {  
public static boolean verbose=false;
 
/*
* prints out the score matrix 
*/
public static void dumpMatrix(int[][] matrix, 
                              String row, 
                              String column){
         
 System.out.print(String.format("%5s",""));
 for (int j =0; j< row.length(); j++){
   System.out.print(String.format("%5s", row.charAt(j)+"  "));
 }          
 System.out.println();
          
for (int i =0; i< column.length(); i++){
  System.out.print(String.format("%5s",column.charAt(i)+ " "));
  for (int j =0; j< row.length(); j++){
     System.out.print(String.format("%5s", matrix[j][i]+" "));
   }
  System.out.println();
 }
}
      
/*
* Needleman-Wunsch Dynamic Programming Algorithm
* see http://amrita.vlab.co.in/?sub=3&brch=274&sim=1431&cnt=1
* Runtime complexity: O(NxM), 
* Space complexity: O(NxM) where N,M are lengths of the sequences
*/
public static AlignmentResult computeNWAlignment(String seq1, 
                                                 String seq2, 
                       AlignmentParameters parameters) {
             
//It is easier to index and faster to use
//integer matrix for fixed length storage
int[][] scoreMatrix = new int[seq1.length()+1][seq2.length()+1];
             
int gapPenalty=parameters.getGapPenalty();
int substitutePenalty=parameters.getSubstitutePenalty();
             
AlignmentResult result = new AlignmentResult();
result.setParameters(parameters);
             
//Initialize the score matrix
//the first row and column are for the gap
//Complexity: O(NxM)
for (int i =0; i< seq1.length()+1; i++){
 for (int j =0; j< seq2.length()+1; j++){
     scoreMatrix[i][j]=0;
     if (i==0) 
    scoreMatrix[i][j] = gapPenalty*j;
     else if (j==0) 
        scoreMatrix[i][j] = gapPenalty*i;
  }
}
 
if (verbose){
   System.out.println("Initial Matrix");
   dumpMatrix(scoreMatrix, "-"+seq1, "-"+seq2);
}
             
int similarityCost=0;
int matchCost=0;
int seq1GapCost=0;
int seq2GapCost=0;
             
//Compute the minimum cost scores between all 
//possible pairs of prefixes
//Complexity: O(NxM)
for (int i =1; i< seq1.length()+1; i++){
for (int j =1; j< seq2.length()+1; j++){
                     
//Case 1: The cost of mistmatch between the two prefixes
similarityCost= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : 
                                   substitutePenalty;   
 
matchCost = scoreMatrix[i-1][j-1] + similarityCost;
                     
//Case 2: the cost of adding a gap on sequence 2
seq2GapCost = scoreMatrix[i-1][j] + gapPenalty;
                     
//Case 3: the cost of adding a gap on sequence 1
seq1GapCost = scoreMatrix[i][j-1] + gapPenalty;
                     
scoreMatrix[i][j] = Math.min(Math.min(matchCost,seq1GapCost),
                             seq2GapCost);
 }
}
             
if (verbose){
 System.out.println("\nFilled Matrix");
 dumpMatrix(scoreMatrix, "-"+seq1, "-"+seq2);
}
             
//Reconstruct the Alignment by backtracking on 
//the score matrix
//Complexity O(N+M)
StringBuilder alignedSequence1= new StringBuilder();
StringBuilder alignedSequence2= new StringBuilder();
             
int j = seq2.length();
int i = seq1.length();
             
while (i >0 || j > 0) {
 if (i>0 && j > 0) 
    similarityCost= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : 
                     substitutePenalty;
    else
    similarityCost = Integer.MAX_VALUE;
                  
if ( i > 0 && j >0 && 
     scoreMatrix[i][j] == scoreMatrix[i-1][j-1] + similarityCost) { 
    alignedSequence1.append(seq1.charAt(i-1));
        alignedSequence2.append(seq2.charAt(j-1));
        i=i-1;
    j=j-1;
 }
 else if ( i > 0 && 
      scoreMatrix[i][j] == scoreMatrix[i-1][j] + gapPenalty){
    alignedSequence2.append("-");
    alignedSequence1.append(seq1.charAt(i-1));
    i=i-1;
 }
else if ( j > 0 && 
       scoreMatrix[i][j] == scoreMatrix[i][j-1] + gapPenalty){
    alignedSequence1.append("-");
    alignedSequence2.append(seq2.charAt(j-1));
    j=j-1;
}
} // end of while
             
 result.setTotalCost(scoreMatrix[seq1.length()][seq2.length()]);
 result.setAlignmentLength(alignedSequence1.length());
 result.setAlignments(new String[] 
    {alignedSequence1.reverse().toString(), 
     alignedSequence2.reverse().toString()});
                 
 return result;
 }

public static double score(AlignmentResult result){
	String[] alignments= result.getAlignments();
	//result
	int matches=0;
	int gaps=0;
//	StringBuffer sb = new StringBuffer();
	for (int k=0; k < alignments[0].length(); k++){
	if (alignments[0].charAt(k)==alignments[1].charAt(k)) {
	    matches++;
	//    sb.append("|");
	} //else sb.append(" ");
	               
	  if ( (alignments[0].charAt(k)=='-') ||
	       (alignments[1].charAt(k)=='-')  ) 
	  gaps++;
	}
	//String sc=String.format( "%5.3g", (float)matches/alignments[0].length()).trim();
	//pw.println(nme+" match_score="+matches+" identity="+sc+" gaps="+gaps+" edit_distance="+result.getTotalCost()+" length="+result.getAlignmentLength());
	//pw.println(alignments[0]);
	//pw.println(sb.toString());
	//pw.println(alignments[1]);
	return (double)matches/(double) alignments[0].length();
}
public static void printResult(PrintWriter pw, AlignmentResult result, String nme){
	String[] alignments= result.getAlignments();
    
	
	int matches=0;
	int gaps=0;
	StringBuffer sb = new StringBuffer();
	for (int k=0; k < alignments[0].length(); k++){
	if (alignments[0].charAt(k)==alignments[1].charAt(k)) {
	    matches++;
	    sb.append("|");
	} else sb.append(" ");
	               
	  if ( (alignments[0].charAt(k)=='-') ||
	       (alignments[1].charAt(k)=='-')  ) 
	  gaps++;
	}
	String sc=String.format( "%5.3g", (float)matches/alignments[0].length()).trim();
	pw.println(nme+" match_score="+matches+" identity="+sc+" gaps="+gaps+" edit_distance="+result.getTotalCost()+" length="+result.getAlignmentLength());
	pw.println(alignments[0]);
	pw.println(sb.toString());
	pw.println(alignments[1]);
}
}