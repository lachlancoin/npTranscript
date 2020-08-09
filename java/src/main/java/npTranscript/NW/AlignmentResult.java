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
     
}