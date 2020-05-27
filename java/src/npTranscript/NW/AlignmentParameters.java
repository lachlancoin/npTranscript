package npTranscript.NW;

/*
 *  Alignment parameters used in Needleman-Wunsch Algorithm 
* with gap penalty
 */
public class AlignmentParameters {
 
    private int gapPenalty;
     
    private int substitutePenalty;
 
    public int getGapPenalty() {
        return gapPenalty;
    }
 
    public void setGapPenalty(int gapPenalty) {
        this.gapPenalty = gapPenalty;
    }
 
    public int getSubstitutePenalty() {
        return substitutePenalty;
    }
 
    public void setSubstitutePenalty(int substitutePenalty) {
        this.substitutePenalty = substitutePenalty;
    }
 
    public AlignmentParameters(int gapPenalty, 
                                  int substitutePenalty) {
        super();
        this.gapPenalty = gapPenalty;
        this.substitutePenalty = substitutePenalty;
    }
     
    public AlignmentParameters() {
        super();
        this.gapPenalty = 2;
        this.substitutePenalty = 1;
    }
}