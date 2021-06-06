package npTranscript.cluster;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

public class EmptyAnnotation extends Annotation {
String chrom1;
	public EmptyAnnotation(String chrom, String string, int seqlen, PrintWriter pw) throws IOException {
		//pw.println("ID\tName\tdescription");
		this.chrom1 = chrom;
		pw.println(chrom1+"\t"+string);
		// TODO Auto-generated constructor stub
	}
	@Override
	public String nextDownstream(int rightBreak,  boolean forward){
		return chrom+"_"+TranscriptUtils.round(rightBreak,CigarHash2.round);
	}
	@Override
	public String  nextUpstream(int rightBreak, boolean forward){
		return chrom+"_"+TranscriptUtils.round(rightBreak,CigarHash2.round);
	}
	public void adjust3UTR(int seqlen2) {
		
	}
	public boolean isLeader(int prev_pos) {
		return prev_pos<100;
	}
	
	

}
