package npTranscript.cluster;

import java.io.File;
import java.io.IOException;

public class EmptyAnnotation extends Annotation {

	public EmptyAnnotation(String chrom, int seqlen) throws IOException {
		super(chrom, seqlen);
		// TODO Auto-generated constructor stub
	}
	public void adjust3UTR(int seqlen2) {
		
	}
	public boolean isLeader(int prev_pos) {
		return false;
	}
	
	public String nextDownstream(int rightBreak, int chrom_index){
		
		return chrom_index+"."+TranscriptUtils.round(rightBreak, CigarHash2.round);
	}
	
	public String nextUpstream(int leftBreak, int chrom_index){
		
		return chrom_index+"."+TranscriptUtils.round(leftBreak, CigarHash2.round);
	}

}
