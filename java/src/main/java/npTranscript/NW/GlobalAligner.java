package npTranscript.NW;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;

import japsa.bio.misc.dnaPlatform.gui.proteinConvertUtilities;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;

public class GlobalAligner {
 public static void main(String[] args){
	 try{
		 String ref="/home/lachlan/github/npTranscript/data/SARS-Cov2/wuhan/wuhan_coronavirus.fasta.gz";
		Sequence reference = SequenceReader.readAll(ref, Alphabet.DNA()).get(0);
		
		 AlignmentParameters ap =new AlignmentParameters();
		 
		 PrintWriter pw ;
		
			 OutputStream os = new FileOutputStream("align.txt");
			 pw =  new PrintWriter(os);
		 String[] args1 = args[0].split(":");
		 String[] args2 = args[1].split(":");
		ArrayList<Sequence>[] seqs = new ArrayList[args1.length];
		for(int i=0; i<seqs.length; i++){
			seqs[i] = SequenceReader.readAll(args1[i], Alphabet.DNA());
		}
		for(int i=0; i<seqs[0].size(); i++){
			Sequence seqs_0 = seqs[0].get(i);
			String nme = seqs_0.getName();
			String[] pos = nme.split("\\.")[1].split("_");
			Sequence subseq = reference.subSequence(Integer.parseInt(pos[0])-1, Integer.parseInt(pos[1]));
			String ref1 = subseq.toString();
			
			CodonTable ct = CodonTableFactory.createUniversalTranslator();
			char[] prot =  translate(ref1.toCharArray(), ct);
	//	System.err.println(">"+nme);
	//	System.err.println(new String(prot));
	//	System.err.println(ref1);
	//	System.err.println("");
			
			Sequence seqs_samp1 = null;
			char[] prot_samp1 = null;
			
			for(int j=0 ;j<seqs.length; j++){
			  Sequence seqs_1 = get(seqs[j], nme);	
			  AlignmentResult result= NeedlemanWunsch.computeNWAlignment(ref1, seqs_1.toString(), ap);
			  String[] alignments = result.getAlignments();
			  char[] res2 = removeGapsInSecond(alignments[1].toCharArray(), alignments[0].toCharArray());
			  char[] prot1 = translate(res2,ct);
			  
			  
			  try{
			  NeedlemanWunsch.printResult(pw, result,expand(prot), expand(prot1),  seqs_1.getName()+": Ref vs "+args2[j]);
			  }catch(Exception exc){
				  exc.printStackTrace();
			  }
			  if(j==0){
				  prot_samp1 = prot1;
				  seqs_samp1 = seqs_1;
			  }else{
				  AlignmentResult result1= NeedlemanWunsch.computeNWAlignment(seqs_samp1.toString(), seqs_1.toString(), ap);
				  NeedlemanWunsch.printResult(pw, result1,expand(prot_samp1), expand(prot1),  seqs_samp1.getName()+": "+args2[0]+" vs "+args2[j]);

			  }
			}
			
			
			 pw.println();
		}
		pw.flush();
		pw.close();
	 }catch(Exception exc ){
		 exc.printStackTrace();
	 }
 }
private static char[] removeGapsInSecond(char[] string, char[] string2) {
	StringBuffer sb = new StringBuffer();
	for(int i=0; i<string.length; i++){
		if(string2[i]!='-'){
			sb.append(string[i]);
		}
	}
	return sb.toString().toCharArray();
}
static String expand(char[] c){
	if(c==null) return null;
//	char[] c = st.toCharArray();
	int len = c.length;
	char[] out = new char[len*3];
	for(int i=0; i<len; i++){
		out[i*3] = '.';
		out[i*3+1]=c[i];
		out[i*3+2] = '.';
	}
	return new String(out);
}
private static char[] translate(char[] dna, CodonTable ct) {
	int len1 = dna.length;
	int rem = (len1 %3);
	if(rem!=0) return null;
	int len_p = (int)Math.floor((double)len1/3.0);
	char[] p = new char[len_p];
	char[] c1 = new char[3];
	for(int i=0; i<len_p; i++){
		c1[0] = dna[i*3];
		c1[1] = dna[i*3 +1];
		c1[2] = dna[i*3+2];
		p[i] = ct.getAminoAcidChar(c1);
	}
	return p;
}

private static Sequence get(ArrayList<Sequence> l, String name) {
	for(int i=0; i<l.size(); i++){
		Sequence seq = l.get(i);
		if(seq.getName().equals(name)){
			return seq;
		}
		
	}
	return null;
}
}
