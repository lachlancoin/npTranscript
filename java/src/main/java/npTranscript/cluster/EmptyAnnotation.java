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
	public EmptyAnnotation( ) throws IOException {
		
		// TODO Auto-generated constructor stub
	}
	@Override
	public int nextDownstream(int rightBreak,  boolean forward, int round){
		return TranscriptUtils.round(rightBreak,round);
	}
	@Override
	public int  nextUpstream(int rightBreak, boolean forward, int round){
		return TranscriptUtils.round(rightBreak,round);
	}
	public void adjust3UTR(int seqlen2) {
		
	}
	public boolean isLeader(int prev_pos) {
		return prev_pos<100;
	}
	
	
	/*if(forward!=null){
		if(strand.charAt(jjk)=='+') secondKey.append(TranscriptUtils.round(breaks.get(e1-1), CigarHash2.round1));
		else secondKey.append(TranscriptUtils.round(breaks.get(s1),CigarHash2.round1));
	}else{
		secondKey.append(TranscriptUtils.round(breaks.get(e1-1), CigarHash2.round1));
	}*/
	
	
	public String nextDownstream(int rightBreak, int chrom_index, Boolean forward, int round){
		return ""+TranscriptUtils.round(rightBreak,round);
	}
	public String  nextUpstream(int rightBreak, Boolean forward, int round){
	//	secondKey.append(this.chrom_);		secondKey.append("/");
	//	secondKey.append(this.strand);secondKey.append("/");
		
		return ""+TranscriptUtils.round(rightBreak,round);
	}

}
