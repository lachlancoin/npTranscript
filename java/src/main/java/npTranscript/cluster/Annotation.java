package npTranscript.cluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.stream.Collectors;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import japsa.seq.Sequence;
import japsa.seq.ZipGFF;
import npTranscript.cluster.Annotation.Interval;
import npTranscript.run.ViralTranscriptAnalysisCmd2;
/**
 * @author Lachlan Coin
 *
 */

public class Annotation{
	int min_overlap=5;

	protected static class Interval{
		int st; int end;
		int len1; int len2; // parent feature lengths
		String gene;
		int length(){
			return end-st;
		}
		public Interval(int st2, int end2, String gene) {
			this(st2,end2, gene,end2-st2, end2-st2);
		}
		
		
		public Interval(int st2, int end2, String gene, int len1, int len2) {
			this.st = st2;
			this.end = end2;
			this.gene = gene;
			this.len1 = len1;
			this.len2 = len2;
		}
		double prop1(){
			return length()/len1;
		}
		double prop2(){
			return length()/len2;
		}
		
		public String toString(){
			double prop1 = (double) length()/(double) len1;
			double prop2 = (double) length()/(double) len2;
			return st+","+end+","+gene+","+String.format("%5.3g", prop1)+","+String.format("%5.3g", prop2);
		}
		public double overlap(int st1, int end1){
			double k0= Math.min(end1-st1, end-st);
			double k1 = Math.min(end1-st, end -st1);
			return Math.min(k0, k1);
		}
		public Interval overlapI(int st1, int end1){
			List<Integer>d = Arrays.asList(new Integer[] {end1-st1,end-st, end1-st, end-st1});
			int min_i = d.indexOf(Collections.min(d));
			int len1 = this.length();
			int len2 = end1 = st1;
			if(min_i==0) return new Interval(st1,end1,gene, len1, len2);
			if(min_i==1) return new Interval(st,end, gene, len1, len2);
			if(min_i==2) return new Interval(st, end1, gene, len1, len2);
			if(min_i==3) return new Interval(st1,end, gene, len1, len2);
			return null;
		}
	}
	protected List<Interval> intervalsF = new ArrayList<Interval>();
	protected List<Interval> intervalsR = new ArrayList<Interval>();
	
		protected List<String> genes= new ArrayList<String>();
		protected  List<Integer> start = new ArrayList<Integer>();
		protected  List<Integer> end = new ArrayList<Integer>();
		protected  List<Boolean> strand = new ArrayList<Boolean>();
		 
		
		 public void clear(){
				this.genes.clear();
				this.start.clear();
				this.end.clear();
				this.strand.clear();
			}
		 
	//	List<Integer> breakSt =  new ArrayList<Integer>();
		//List<Integer> breakEnd =  new ArrayList<Integer>();
		//  int seqlen;

		//Name,Type,Minimum,Maximum,Length,Direction,gene
	//	5'UTR,5'UTR,1,265,265,forward,none

		public static String NAstring = "NA";
		public static int tolerance = 10;
		public static int correctionDistLeft = 10000;
		public static int correctionDistRight = 10000;	
		public static boolean enforceStrand = true;
		 protected int seqlen;
		public void print(PrintWriter pw){
			
			int len = spliced_count.length;
			pw.print("Gene\tStart\tEnd");
			for(int j=0; j<len; j++){
				pw.print("\tSpliced_5_"+j);
				pw.print("\tUnspliced_5_"+j);
			}
			pw.println();
			for(int i=0; i<start.size(); i++){
				pw.print(genes.get(i)+"\t"+start.get(i)+"\t"+end.get(i));
				for(int j=0; j<len; j++){
					pw.print("\t"+spliced_count[j][i]);
					pw.print("\t"+unspliced_count[j][i]);
				}
				pw.println();
			}
			pw.flush();
				}
				static String endStr="end";
		public String  nextDownstreamG(int rightBreak,  boolean forward){
			int i=nextDownstream(rightBreak, forward);
			if(i==genes.size()) {
				return genes.get(i-1);
			}
			return this.genes.get(i);
		}
			
		
		public int nextDownstream(int rightBreak,  boolean forward){
		
			for(int i=0; i<start.size(); i++){
				if(!enforceStrand || forward==this.strand.get(i)){
					if(rightBreak -tolerance <= start.get(i) ){//&& rightBreak < end.get(i)){
						return i;
					}
				}
			}
			return genes.size();//+chrom;//+(chrom_index>0 ? "."+chrom_index : "");//+ (enforceStrand ? (forward ? "+" : "-") : "");
		}
		/*public String convert(List<Integer> l, int chrom_index){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<l.size(); i++){
				if(l.get(i)==null) sb.append("null");
				else if(l.get(i)<0) sb.append("end"+(chrom_index>0 ? "."+chrom_index : ""));//+ (enforceStrand ? (forward ? "+" : "-") : "");)
				else if(l.get(i)>=start.size()) sb.append("start"+(chrom_index>0 ? "."+chrom_index : ""));//+ (enforceStrand ? (forward ? "+" : "-") : "");)
				else sb.append(genes.get(i));
				if(i%2 ==0){
					sb.append("_") ;
				}else sb.append(","); 
			}
			return sb.toString();
		}*/
		
		public String nextUpstreamG(int leftBreak, boolean forward){
			int i = this.nextUpstream(leftBreak, forward);
			return genes.get(i);
		}
		public int nextUpstream(int leftBreak,  boolean forward){
		//	if(leftBreak<0) return  "null";
			for(int i=start.size()-1; i>=0 ;i--){
				if(!enforceStrand || forward==this.strand.get(i)){
				if(leftBreak+tolerance >= start.get(i)){// && leftBreak-tolerance<end.get(i)){
					return i;
				}
				}
			}
			return -1;//+chrom;//+(chrom_index>0 ? "."+chrom_index : "");//+(enforceStrand ? (forward ? "+" : "-") : "");
		}
	
    int updateLeft(int en, int refStart){
    	int breaks_out_0= en;
    	for(int i=0; i< this.start.size(); i++){
    		int st_i = start.get(i);
	    	if(en+ tolerance > st_i ){
				if(Math.abs(st_i - en) < correctionDistLeft && st_i >refStart) {
					breaks_out_0 = st_i;
					
				}
			}
    	}
    	return breaks_out_0;
    }
    
    int updateRight( int st, int refEnd){
    	for(int i=0; i< this.start.size(); i++){
	    	int st_i = start.get(i);
	    	if(st - tolerance < st_i) {
				int dist = Math.abs(st_i - st);
				if(dist < correctionDistRight && st_i <refEnd) {
					return st_i;
				}else{
					return st;
				}
			}
    	}
    	return st;
    }
		
	public void getBreaks(int[] breaks_in, int[] breaks_out, int refStart, int refEnd){
		System.arraycopy(breaks_in, 0, breaks_out, 0, 2);
		breaks_out[1] = this.updateRight(breaks_in[1], refEnd);
		breaks_out[0] = this.updateLeft(breaks_in[0], refStart);
	}
	public Annotation(){
		
	}
	

	public void updateChrom( String chrname) {
		this.chrom = chrname;
	}
	
	public static Annotation getAnnotation(String annot_file,  PrintWriter annotation_pw, int source_count) throws ZipException, IOException{
		ZipFile  anno  = null;
		File gffFile = new File(annot_file);
		if(!gffFile.exists()) gffFile = null;
		boolean writeDirect = true;
		if(gffFile!=null && (gffFile.getName().indexOf(".gff")>=0 || gffFile.getName().indexOf(".gtf")>=0)){
			if(gffFile.getName().endsWith(".zip")){
				anno = new ZipFile(gffFile);
				System.err.println(anno.getName());
			}
			else {
				String out_nme = gffFile.getName();
				int ind = out_nme.lastIndexOf('.');
				out_nme = out_nme.substring(0, ind);
				File outzip = new File(gffFile.getParentFile(),out_nme+".zip");
				if(outzip.exists()){
					System.err.println("reading existing gff zip file "+outzip.getAbsolutePath());
					anno = new ZipFile(outzip);
				}
				else{
					System.err.println("making gff.zip file");
					ZipGFF gffin =  new ZipGFF( gffFile, outzip,writeDirect);
					gffin.run();
					anno = new ZipFile(outzip);
				}
			}
		}
		
		
		Annotation annot  = null;
		if(gffFile!=null && (gffFile.getName().indexOf(".gff")>=0 || gffFile.getName().indexOf(".gtf")>=0)){
				annot = new GFFAnnotation(gffFile,annotation_pw, gffFile.getName().indexOf(".gff")<0);
				
		}else{
			
			annot=	new Annotation(new File(annot_file),  annotation_pw, source_count);
		}
		return annot;
	
		
	
	}
	
	public void close(){
		
	}
	
	public	Annotation(File f, int source_count) throws IOException{
		this(f, null, source_count);
	}
	
	public	Annotation(File f,PrintWriter pw, int source_count) throws IOException{
		//this(chrom, seqlen);
		BufferedReader br = new BufferedReader(new FileReader(f)) ;
			List<String> header = Arrays.asList(br.readLine().split(","));
			int gene_ind = header.indexOf("gene");
			int start_ind = header.indexOf("Minimum");
			int end_ind = header.indexOf("Maximum");
			int direction_ind = header.indexOf("Direction");
			String str = "";
			for(int i=0; (str=br.readLine())!=null; i++){
				String[] line = str.split(",");
				boolean both = line[direction_ind].equals("both");
				boolean forward = line[direction_ind].equals("forward");
				String gene = line[gene_ind];	
				
				int st = Integer.parseInt(line[start_ind]);
				int en = Integer.parseInt(line[end_ind]);
				if(pw!=null){
					String str1 = "0\t"+gene+"\t"+gene+"\t"+gene+"\t"+gene+"\t"+gene;
					pw.println(str1);
				}
				genes.add(gene);
				start.add(st);
				end.add(en);
				Interval interval = new Interval(st,en,gene);
				if(both && enforceStrand){
					strand.add(true); strand.add(false);
					genes.add(gene);
					start.add(st);
					end.add(en);
					intervalsF.add(interval) ;
					intervalsF.add(interval) ;
					
				}else{
					strand.add(forward);
					if(forward) {
						intervalsF.add(interval) ;
						}else {
							intervalsR.add(interval);
						}
				}
			
				//	if(i>0){
				//		this.breakSt.add(break_start);
				//		this.breakEnd.add(st);
				//	}
					//break_start = st;
			
			}
			br.close();
			mkOrfs(source_count);
			this.seqlen = this.end.get(end.size()-1);
		}
	
	String chrom="";
	 protected void mkOrfs(int source_count) {
		orfs = new int[start.size()][];
		overlap = new double[start.size()];
		orf_len = new int[start.size()];
		this.unspliced_count = new int[source_count][start.size()];
		this.spliced_count = new int[source_count][start.size()];
		for(int i=0; i<orf_len.length; i++) {
			orf_len[i] = end.get(i) - start.get(i)+1;
			orfs[i] = new int[] {start.get(i),end.get(i)};
		}
		
	}
		private int[][] orfs; 
		public int[][] unspliced_count;
		public int[][] spliced_count;

	public  double[] overlap;
	public  int[] orf_len;
	
	/*	public int seqlen() {
			return seqlen;
		}*/

		public synchronized void adjust3UTR(int seqlen2) {
			if(end.size()==0) return ;
			if(this.end.get(end.size()-1)>seqlen2 && this.start.get(end.size()-1) < seqlen2){
				System.err.println("adjusting 3UTR end "+seqlen2);
				end.set(end.size()-1, seqlen2);
			}
			
		}

		public boolean isLeader(int prev_pos) {
			return prev_pos < this.start.get(1)+tolerance;
		}

	
		public void addCount(int i, int source_index, boolean spliced) {
			if(spliced) {
				spliced_count[source_index][i]+=1;
				//if(i==2){
				//	System.err.println("hh");
			//	}
			}
			else unspliced_count[source_index][i]+=1;
			
		}
		/*public String getTypeNme(int start, int end, boolean forward) {
			return nmes[getTypeInd(start,end, forward)];
			//return null;
		}*/
		/*public int getTypeInd(int start, int end, Boolean forward) {
			if(start <=TranscriptUtils.startThresh) return end >= seqlen -TranscriptUtils.endThresh ? 0:1;
			else return end >= seqlen -TranscriptUtils.endThresh ? 2:3;
			//return null;
		}*/

		public static String[] nmes = new String[] {"5_3", "5_no3", "no5_3", "no5_no3"};
		public  char getStrand(Iterator<Integer> l){
			if(!l.hasNext()) return 'N';
			Boolean strandi = this.strand.get(l.next());
			while(l.hasNext()){
				if(!strand.get(l.next()).equals(strandi)){
					return 'N';
				}
			}
			return strandi ? '+' : '-';
			
		}

		public String getSpan(List<Integer>breaks, Boolean forward, Collection<Integer> span, SortedSet<String> parents) {
			// TODO Auto-generated method stub
			return ".";
		}

		public String getString(Collection<Integer> span, SortedSet<String> geneNmes) {
			if(span.size()==0) return "-";
			StringBuffer sb = new StringBuffer();
			for(Iterator<Integer> it = span.iterator(); it.hasNext();){
				String gene = this.genes.get(it.next());
				geneNmes.add(gene);
				sb.append(gene);
				if(it.hasNext())sb.append(";");
			}
			return sb.toString();
		}

		
		public void annotate(CigarCluster nxt) {
			/*if(ViralTranscriptAnalysisCmd2.coronavirus){
				if(includeStartEnd || breaks.size()==2){
					secondKey.append(annot.nextUpstream(startPos,chrom_index, forward)+delim);
				}
				boolean firstBreak=true;
				for(int i=1; i<breaks.size()-1; i+=2){
					int gap = breaks.get(i+1)-breaks.get(i);
					String upst = annot.nextUpstream(breaks.get(i), chrom_index,forward); //5prime break
					secondKey.append(upst+delim1);
					secondKey.append(annot.nextDownstream(breaks.get(i+1), chrom_index,forward)+delim);  //3prime break
					if(gap > break_thresh){
						if(firstBreak){
							firstBreak=false;
							if(annotByBreakPosition) parent.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
						//	if(bp!=null) this.bp.addBreakPoint(source_index, 0, breaks.get(i), breaks.get(i+1));
							hasLeaderBreak = true;
						}else{
							if(annotByBreakPosition)  parent.addBreakPoint(source_index, 1, breaks.get(i), breaks.get(i+1));

						}
					}
			
				}
				if(includeStartEnd || breaks.size()==2) {
					secondKey.append(annot.nextUpstream(breaks.get(breaks.size()-1), chrom_index, forward)); //last break is upstream start pos
				}
					
			}else */
			throw new RuntimeException ("need to re-implement");/// TODO Auto-generated method stub
			
		}

		protected static List<String> empty_list = Arrays.asList(new String[0]);

		protected boolean matchChrom(String chrom) {
			// TODO Auto-generated method stub
			return this.chrom.equals(chrom);
		}
		public List<String> matchExons(List<Integer> br, String chrom, Boolean forward) {
			if(br.get(0).equals(612)){
				if(br.get(1).equals(1884)){
					System.err.println("h");
				}
			}
			List<Integer> l = new ArrayList<Integer>();
			if(!matchChrom(chrom)){
				return empty_list;
			}
			List<Interval> intervals =( !enforceStrand || forward)? intervalsF : intervalsR;
List<String> genes = new ArrayList<String>();
			for(int i=0;i<br.size(); i+=2){
				int st1 = br.get(i);
				int end1 = br.get(i+1);

				List<Interval> overlaps=intervals.stream().filter(c -> c.overlap(st1, end1) > min_overlap).collect(Collectors.toList());
				List<Interval> overlap1 = overlaps.stream().map(c->c.overlapI(st1, end1)).collect(Collectors.toList());
				List<String> genes1=	overlaps.stream().map(c->c.gene).collect(Collectors.toList());
                if(genes1.size()==0){
                	genes.add("NA");genes.add("NA");
                }else if(genes1.size()==1){
                	genes.add(genes1.get(0));genes.add(genes1.get(0));
                }else{
                	genes.add(genes1.get(0));
                	genes.add(genes1.get(genes1.size()-1));
                }
			}
			
//			List<String> genes=	intervals.stream().filter(c -> c.overlap(st1, end1) > min_overlap).map(c->c.gene).collect(Collectors.toList());
			return genes;
		}

		
		public List<String> matchExons1(List<Integer> br, String chrom, Boolean forward) {
			List<Integer> l = new ArrayList<Integer>();
			if(!this.chrom.equals(chrom)){
				return empty_list;
			}
			l.add(this.nextUpstream(br.get(0) , forward));
			for(int i=1;i<br.size()-1; i++){
					if(i%2 ==0 ){//exon start
						l.add(this.nextDownstream(br.get(i) , forward));
					}else{ //exon end
						l.add(this.nextUpstream(br.get(i) , forward));
					}
				}
			l.add(this.nextDownstream(br.get(br.size()-1) , forward));

			
			return getStr(l);
		}



		protected List<String> getStr(List<Integer> l) {
			String[] res = new String[l.size()];
			for(int i=0; i<res.length; i++)res[i] = l.get(i) <0 ? "start" : l.get(i)==genes.size() ?  "end":  this.genes.get(l.get(i));
			return Arrays.asList(res);
		}


		public int seqLen() {
			return seqlen;
		}


		
	}