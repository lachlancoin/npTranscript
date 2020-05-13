package npTranscript.run;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import japsa.bio.np.barcode.SWGAlignment;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.CigarHash2;

@Deployable(scriptName = "npConsensus.run", scriptDesc = "Analysis of consensus")
public class ConsensusMapper extends CommandLine {

	static SequenceOutputStream os;
	
	public ConsensusMapper() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("fasta", null, "Name of consensus fasta file", true);
		addString("reference", null, "Name of reference genome", true);
		addInt("offset",0, "offset to add from breakpoints",false); // this included for earlier versions of npTranscript
		addStdHelp();

	}
	public static void main(String[] args1) throws IOException, InterruptedException {
		CommandLine cmdLine = new ConsensusMapper();
		String[] args = cmdLine.stdParseLine(args1);

		String fasta = cmdLine.getStringVal("fasta");
		String reference = cmdLine.getStringVal("reference");
		int offset = cmdLine.getIntVal("offset");
		run(fasta, reference, offset);
	}
	static Pattern p1 = Pattern.compile("-{1,}"); 
	static Pattern p2 = Pattern.compile("-{2,}"); //consider two or more to be deletion

	static Pattern p3 = Pattern.compile("-{3,}"); //consider two or more to be deletion
	
	static class Indel{
		int start;
		String seq;
		//int len;
		Indel(int start, String seq){
			this.start = start;
			this.seq = seq;
			//this.len = len;
		}
		public String toString(){
			//if(len==1) return start+"";
			//else
				return start+":"+seq;
		}
		public int len() {
			// TODO Auto-generated method stub
			return seq.length();
		}
	}
	//getString(List<Indel>)
	//gaps in seq ; non gaps in seq2
	static int getGaps(char[] seq,char[] seq2, Pattern p, List<Indel> allMatches){
		String seqr = new String(seq);
		String seqr2 = new String(seq2);
		Matcher m = p.matcher(seqr);
		int tot =0;
		while (m.find()) {
				int len =  m.end()-m.start();
				tot+=len;
			   allMatches.add(new Indel(m.start(),seqr2.substring(m.start(), m.end())));
		}
		return tot;
	}
	
	/** removes positions in first sequence which cause gaps in second sequence */
	private static String removeGapsInRef(char[] seq1, char[]  seqRef,List<Indel> gaps, int startInRef) {
		if(seq1.length!=seqRef.length) throw new RuntimeException("!!");
		gaps.clear();
		int gap = getGaps(seqRef,seq1,p1, gaps);
		char[] seq1_new = new char[seq1.length - gap];
		int st =0;
		int target_st =0;
		int total_gap =0;
		for(int i=0; i<gaps.size(); i++){
			Indel indel = gaps.get(i);
			int end = indel.start; // end of previous non gapped sequence
			int len = end-st;
			if(len>0){
				for(int j=0; j<len; j++){
					seq1_new[target_st +j] = seq1[st+j]==seqRef[st+j] ? seq1[st+j] : lower(seq1[st+j]);
				}
				//System.arraycopy(seq1, st, seq1_new, target_st, len);
				target_st+=len;
			}
			st = end+ indel.len();
			indel.start = indel.start-total_gap+startInRef;  // adjust indel so that start position now is in ref coords (i.e. remove total_gap
			total_gap+=indel.len();
		}
		if(st < seqRef.length){
			int len = seqRef.length - st;
			for(int j=0; j<len; j++){
				
				seq1_new[target_st +j] = seq1[st+j]==seqRef[st+j] ? seq1[st+j] : lower(seq1[st+j]);
			}
			//System.arraycopy(seq1, st, seq1_new, target_st, len);
		}
		return new String(seq1_new);
	}
	
	static List<Integer> readBreaks(String[]br){
		Integer[] breaks = new Integer[br.length];
		for(int k=0; k<br.length; k++) breaks[k] = Integer.parseInt(br[k]);
		return Arrays.asList(breaks);
	}
	
	private static char lower(char c) {
		return (new String(new char[] {c})).toLowerCase().charAt(0);
	}
	private static void run(String fasta, String refFile, int offset) throws IOException {
		// TODO Auto-generated method stub
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());
		List<String> header = Arrays.asList(
				"ID      chrom start   end     type_nme        startBreak      endBreak        isoforms        leftGene        rightGene       totLen  countTotal".split("\\s+")
		);
		int genome_index = Integer.parseInt(fasta.split("\\.")[0]);
		Sequence refSeq = genomes.get(genome_index);
		ArrayList<Sequence> reads = SequenceReader.readAll(fasta, Alphabet.DNA());
		
		File fasta_in = new File(fasta);
		String fastanme = fasta_in.getName().replaceAll(".fasta", "").replaceAll(".fa", ".fa");
		File fasta_out = new File(fasta_in.getParentFile(), fastanme+".ref_aligned.fa");
		OutputStreamWriter os1 = new OutputStreamWriter(new FileOutputStream(fasta_out));
		
		File fasta_out1 = new File(fasta_in.getParentFile(), fastanme+".ref.fa");
		OutputStreamWriter os2 = new OutputStreamWriter(new FileOutputStream(fasta_out1));
		
		
		File fasta_out2 = new File(fasta_in.getParentFile(),fastanme+".unmapped_insertion.fa");
		OutputStreamWriter os3 = new OutputStreamWriter(new FileOutputStream(fasta_out2));
	
		//SequenceOutputStream os =new SequenceOutputStream(os1);
		int step = 100;
		for(int i=0; i< reads.size(); i++){
			Sequence readSeq = reads.get(i);
			String nme = readSeq.getName();
			String[] desc = readSeq.getDesc().split(";"); //ID0.0;MT007544.1;1;29893;5_3;22;27904;2;leader;ORF8;-1;2262
			List<Integer> breaks = ConsensusMapper.readBreaks(desc[0].split(","));
			
	//		int startBr = Integer.parseInt(desc[0].split(",")[0]); 
	//		int endBr = Integer.parseInt(desc[0].split(",")[1]); 
			//int startBr = Integer.parseInt(desc[header.indexOf("startBreak")]);
			//int endBr = Integer.parseInt(desc[header.indexOf("endBreak")]);
			//String leftGene = 			desc[header.indexOf("leftGene")];
			//String rightGene = 			desc[header.indexOf("rightGene")];
			
			
			
			int[] start1 = new int[(int) Math.floor((double)breaks.size()/2.0)];
			int[] end1= new int[start1.length];
			int[] start = new int[start1.length];
			int[] end = new int[ end1.length];
			//double[] ident = new double[start1.length];
			List<Indel> [] gaps = new List[start1.length];
			String[] sequ = new String[start1.length];
			String[] remaining = new String[start1.length];
			String[] sequ_refonly = new String[start1.length];
			SWGAlignment[] align = new SWGAlignment[start1.length];
			List<Indel> allgaps = new ArrayList<Indel>();
			Integer[] breaks_reads = new Integer[breaks.size()];
			Integer[] breaks_ref = new Integer[breaks.size()];
			double ident = 0;
			double totlen = 0;
			for(int k1=start1.length-1; k1>=0; k1--){
				Sequence readSeq1 =k1==start1.length-1 ? readSeq: readSeq.subSequence(0, start1[k1+1]-1);
				int endBr = Math.min( breaks.get(2*k1+1) + offset, refSeq.length());
				int stBr = Math.max(0,breaks.get(2*k1) - offset);
				int right_offset = endBr;
				align[k1] = SWGAlignment.align(readSeq1, refSeq.subSequence(stBr, endBr));
				start1[k1] = align[k1].getStart1();
				start[k1] = align[k1].getStart2();
				end1[k1] = start1[k1] + align[k1].getSequence1().length - align[k1].getGaps1();
				end[k1] = start[k1] + align[k1].getSequence2().length - align[k1].getGaps2();
				
				gaps[k1] = new ArrayList<Indel>();
				sequ[k1] = removeGapsInRef(align[k1].getSequence1(), align[k1].getSequence2(),gaps[k1], start[k1]);
				allgaps.addAll(gaps[k1]);
				sequ_refonly[k1]= (new String(align[k1].getSequence2())).replaceAll("-", "");
				if(sequ[k1].length()!=sequ_refonly[k1].length()) throw new RuntimeException("this should not happen");
				breaks_ref[2*k1] = start[k1]+stBr;
				breaks_ref[2*k1+1] =end[k1]+stBr;
				breaks_reads[2*k1] = start1[k1];
				breaks_reads[2*k1+1] = end1[k1];
				ident+= align[k1].getIdentity();
				totlen +=sequ_refonly[k1].length();
				if(k1<start1.length-1){
					if(start1[k1+1] <= end1[k1]) {
						throw new RuntimeException("problem");
					}else{
						if(end1[k1]+1<start1[k1+1]){
							remaining[k1]=new String(readSeq.subSequence(end1[k1], start1[k1+1]).charSequence());
						}
					}
				}else{
				/*	if(end1[k1]+1<readSeq.length()){
						remaining[k1]=new String(readSeq.subSequence(end1[k1], readSeq.length()).charSequence());
					}*/
				}
			
			}
			String identity =String.format("%5.3g", ident/totlen);
		  writeFasta(os1, sequ, readSeq.getName(),getString(breaks_ref)+";"+getString(breaks_reads)+";"+identity
				 +" "+allgaps.toString().replaceAll("\\s+", ""));
		  writeFasta(os2, sequ_refonly,readSeq.getName(),getString(breaks_ref)+";"+getString(breaks_reads)+";"+identity);
		  writeFasta(os3, remaining,readSeq.getName(),getString(breaks_ref)+";"+getString(breaks_reads));
		/*  if(right_start1 > left_end1+1){
				// remaining = null;
				 Sequence remaining=readSeqLeft.subSequence(left_end1, readSeqLeft.length());
				  writeFasta(os3, new String(remaining.charSequence()),null,  readSeq.getName(),left_start+","+left_end+","+right_start+","+right_end+" "
							+left_start1+","+left_end1+","+right_start1+","+right_end1
									  +" "+readSeq.getDesc());
				//System.err.println("unmapped insertion: "+left_start1+","+left_end1+","+right_start1+","+right_end1+""+remaining);
			}*/
		}
		os1.close();
		os2.close();
		os3.close();
	}
	private static void writeFasta(OutputStreamWriter os1, String[] sequ1,  String name, String desc) throws IOException{
		boolean allnull = true;
		for(int i=0; i<sequ1.length; i++){
			if(sequ1[i]!=null) allnull=false;
		}
		if(allnull) return;
		os1.write(">"+name+" "+desc+"\n");
		int step = 100;
		for(int k=0; k<sequ1.length; k++){
			if(sequ1[k]!=null){
			for(int j=0; j<sequ1[k].length(); j+=step){
				String substr = sequ1[k].substring(j, Math.min(sequ1[k].length(), j+step));
				os1.write(substr);
				os1.write("\n");
			}
			}
		}
	}
	public static String getString(Integer[] l){
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<l.length; i++){
			if(i>0) sb.append(",");
			sb.append(l[i]);
		}
		return sb.toString();
	}
	
}
