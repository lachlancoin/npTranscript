package npTranscript.NW;

/*
 * Class that encapsulates results of a binary Alignment: 
* parameters, edit distance and the actual alignments
 */
public class AlignmentResult {
     
    private int totalCost=0;
     
    private int alignmentLength=0;
     
     
    public int getAlignmentLength() {
        return alignmentLength;
    }
 
    public void setAlignmentLength(int alignmentLength) {
        this.alignmentLength = alignmentLength;
    }
 
     
    private int matches=0;
     
    private AlignmentParameters parameters=null;
     
    private String[] alignments=null;
 
    public int getMatches() {
        return matches;
    }
 
    public void setMatches(int matches) {
        this.matches = matches;
    }
 
    public int getTotalCost() {
        return totalCost;
    }
 
    public void setTotalCost(int totalCost) {
        this.totalCost = totalCost;
    }
 
    public AlignmentParameters getParameters() {
        return parameters;
    }
 
    public void setParameters(AlignmentParameters parameters) {
        this.parameters = parameters;
    }
 
    public String[] getAlignments() {
        return alignments;
    }
 
    public void setAlignments(String[] alignments) {
        this.alignments = alignments;
    }
     
    private int alignmentScore(String seq1, String seq2){
          int totalCost=0;
           
          for (int k=0; k < seq1.length(); k++){
              if (seq1.charAt(k)!=seq2.charAt(k)) {
                    if ( (seq1.charAt(k)!='-') && 
                                (seq2.charAt(k)!='-')  ) 
                            totalCost++;
              }           
              if ( (seq1.charAt(k)=='-') || 
                               (seq2.charAt(k)=='-')  ) {
                 totalCost+=2;
              }
          }
          return totalCost;
    }

	public double getMatchPerc() {
		int matches=0;
		int gaps=0;
//		StringBuffer sb = new StringBuffer();
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
     
}