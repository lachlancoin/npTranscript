package npTranscript.NW;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

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
		 OutputStream pw_os = new FileOutputStream("align.vcf");
		 PrintWriter pw_vcf = new PrintWriter(new OutputStreamWriter(pw_os));
		
			 OutputStream os = new FileOutputStream("align.txt");
			 pw =  new PrintWriter(os);
		 String[] args1 = args[0].split(":");
		 String[] args2 = args[1].split(":");
		 Map<Integer, Sequence[]> seqmap = new TreeMap<Integer, Sequence[]>();
		 int seqlen = args1.length;
	//	ArrayList<Sequence>[] seqs = new ArrayList[args1.length];
		 pw_vcf.println("#pos,ref,type,identical,"+args[1].replace(":", ",")+",gene_nme,aa_pos,aa_phase,aa_ref,"+args[1].replace(":", ","));
		for(int i=0; i<seqlen; i++){
			List<Sequence> l = SequenceReader.readAll(args1[i], Alphabet.DNA());
			for(int j=0; j<l.size(); j++){
				Sequence seql  = l.get(j);
				String nme = seql.getName();
				String[] pos = nme.split("\\.")[1].split("_");
				int start = Integer.parseInt(pos[0])-1;
				int end = Integer.parseInt(pos[1]);
				Sequence[] seqs = seqmap.get(start);
				if(seqs==null) seqmap.put(start,seqs =  new Sequence[seqlen]);
				seqs[i] = seql;
			}
		}
		for(Iterator<Sequence[]> it_s = seqmap.values().iterator(); it_s.hasNext();){
			
			Sequence[] seqs = it_s.next();
			String nme =seqs[0].getName();

			System.err.println(nme);
			String gene_nme = nme.split("\\.")[0];
			String[] pos = nme.split("\\.")[1].split("_");
			int start = Integer.parseInt(pos[0])-1;
			int end = Integer.parseInt(pos[1]);
			Sequence subseq = reference.subSequence(start, end);
			String ref1 = subseq.toString();
			
			CodonTable ct = CodonTableFactory.createUniversalTranslator();
			char[] prot =  translate(ref1.toCharArray(), ct);
	
			
			Sequence seqs_samp1 = null;
			char[] prot_samp1 = null;
			
			char[][] aligns = new char[seqs.length][];
			char[][] prot_align = new char[seqs.length][];
			SortedMap<Integer, int[]> insertions =  new TreeMap<Integer, int[]>();
			for(int j=0 ;j<seqs.length; j++){
			  Sequence seqs_1 = seqs[j];
			//  System.err.println(seqs_1.length()+" "+ref1.length());
			  AlignmentResult result= NeedlemanWunsch.computeNWAlignment(ref1, seqs_1.toString(), ap);
			  String[] alignments = result.getAlignments();
			  
			//   insertions[j]=
			  char[] res2 = removeGapsInSecond(alignments[1].toCharArray(), alignments[0].toCharArray(), insertions, j, seqs.length);
			  char[] prot1 = translate(res2,ct);
			  aligns[j] = res2;
			  prot_align[j] = prot1;
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
			char[] ref1_seq = ref1.toCharArray();
			  writeVCF(pw_vcf, aligns, ref1_seq, prot_align, prot,  start, insertions, gene_nme);
			
			
			 pw.println();
		}
		pw.flush();
		pw.close();
		pw_vcf.close();
	 }catch(Exception exc ){
		 exc.printStackTrace();
	 }
 }
private static void writeVCF(PrintWriter pw_vcf, char[][] align, char[] ref, char[][] prot_align,  char[] prot_ref,int start, SortedMap<Integer, int[]> insertions, String gene) {
	int len = ref.length;
	int len1 = align.length;
	char[] res = new char[2*len1];
	Arrays.fill(res,',');
	char[] prot_res = new char[2*len1];
	Arrays.fill(prot_res,',');
	if(align.length>1 && align[0].length!=align[1].length) throw new RuntimeException ("!!");
	for(int i=0; i<len; i++){
		boolean allsame=true;
		int rem = i %3;
		int pos_aa =Math.floorDiv(i, 3);// Math.floor((double)i/3.0);
		boolean ident = true;
		for(int j=0; j<len1; j++){
			res[2*j+1] = align[j][i];
			if(prot_align[j]!=null){
			prot_res[2*j+1] = prot_align[j][pos_aa];
			}
			if(align[j][i]!=ref[i]){
				allsame=false;
			}
			if(align[j][i]!=align[0][i]) ident=false;
			
		}
		if(!allsame){//variant position
			int pos = start+i+1;
			String type="mut";
			if(res[1]=='-') type="del";
		//	prot_ref[pos_aa];
			String aa_ref = prot_ref==null ? "": prot_ref[pos_aa]+"";
			pw_vcf.println(pos+","+ref[i]+","+type+","+ident+","+new String(res)+","+gene+","+(1+pos_aa)+","+rem+","+aa_ref+new String(prot_res));
		}
	}
//	int[
	Arrays.fill(prot_res,',');

	for(Iterator<Integer> it = insertions.keySet().iterator(); it.hasNext();){
		Integer i =  it.next();
		int[] val = insertions.get(i);
		int pos = start+i;
		int rem = i %3;
		int pos_aa =Math.floorDiv(i, 3);
		String aa_ref = prot_ref==null ? "": prot_ref[pos_aa]+"";
		boolean ident =true;
		for(int k=0; k<val.length; k++){
			if(val[k]!=val[0]) ident=false;
		}
		pw_vcf.println(pos+","+ref[i]+",ins,"+ident+","+paste(val)+","+gene+","+(1+pos_aa)+","+rem+","+aa_ref+new String(prot_res));

	 }
	
}

private static String paste(int[] val) {
	StringBuffer sb = new StringBuffer();
	for(int i=0 ; i<val.length; i++){
		sb.append(","+val[i]);
	}
	return sb.toString();
}
private static char[] removeGapsInSecond(char[] string, char[] string2, SortedMap<Integer, int[]> insertions, int index, int len) {
	StringBuffer sb = new StringBuffer();
	int pos=0;
	for(int i=0; i<string.length; i++){
		if(string2[i]!='-'){
			sb.append(string[i]);
			pos++;
		}else{
		int[] insp = 	insertions.get(pos);
		if(insp==null)	insertions.put(pos, insp = new int[len]);
		insp[index]++;
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
