/*****************************************************************************
 * Copyright (c)Lachlan Coin University of Melbourne, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/*                                                    
 * 01/03/2020 - Lachlan Coin                                       
 ****************************************************************************/

package npTranscript.run;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Executors;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.Annotation;
import npTranscript.cluster.CigarCluster;
import npTranscript.cluster.CigarHash;
import npTranscript.cluster.CigarHash2;
import npTranscript.cluster.EmptyAnnotation;
import npTranscript.cluster.GFFAnnotation;
import npTranscript.cluster.IdentityProfile1;
import npTranscript.cluster.Outputs;
import npTranscript.cluster.TranscriptUtils;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class ViralTranscriptAnalysisCmd2 extends CommandLine {
	static String src_tag = "SC";
	static String pool_tag = "PT";
private static final class CombinedIterator implements Iterator<SAMRecord> {
		private final SAMRecordIterator[] samIters;
		int current_sam_index =0;
		int currentIndex =0; //relates to chromosome
		private final SAMRecord[] currentVals;
		private final boolean[] returned;
		private final int[] cnts ;
		int max;
		Collection<String>[] readList ;
		Set<Integer> chrs;
		private CombinedIterator(SAMRecordIterator[] samIters, int max, Collection<String>[]readList, Set<Integer> chrs) {
			this.samIters = samIters;
			this.readList = readList;
			currentVals = new SAMRecord[samIters.length];
			returned = new boolean [samIters.length];
			this.max = max;
			cnts = new int[samIters.length];
			this.chrs = chrs;
		}

		@Override
		public boolean hasNext() {
			for(int i=0; i<samIters.length; i++){
				if(samIters[i].hasNext() && cnts[i]<max) return true;
			}
			return false;
		}
		int pool_ind=-1;
		private SAMRecord next(int i){
			SAMRecord sr;
			if(currentVals[i]!=null && !returned[i]){
				sr = currentVals[i];
			}else{
				if(cnts[i]<max){
					
					//sr  = samIters[i].next();
				/*	while(!(
							sr==null || 
							(readList==null || readList.contains(sr.getReadName())) ||
							(chrs==null || chrs.contains(sr.getReferenceIndex()))
							)
							)*/
					inner: while(true)
					{
						sr  = samIters[i].next();
						if(sr==null) break inner;
						pool_ind =-1;
						if(readList!=null){
							inner1: for(int i2=0; i2<readList.length; i2++){
								if(readList[i2].contains(sr.getReadName())){
									pool_ind = i2;
									break inner1;
								}
							}
						}
						if(readList!=null && pool_ind>=0) break inner;
						if(readList==null &&  chrs!=null && chrs.contains(sr.getReferenceIndex())) break inner;
						if(readList==null  && chrs==null) break inner;
					}
					if(sr!=null) {
						sr.setAttribute(pool_tag, pool_ind);
						sr.setAttribute(src_tag, i);
					}
					currentVals[i] = sr;
					returned[i] = false; 
				}else{
					sr = null;
				}
			}
			return sr;
		}

		@Override
		public SAMRecord next() {
			//int curr_index = current_sam_index;
			SAMRecord sr = next(current_sam_index);
			
			if(sr==null || sr.getReferenceIndex()>currentIndex ){
				int[] ref_inds = new int[samIters.length];
				int min_ind =-1;
				int minv = Integer.MAX_VALUE;
				for(int i=0; i<samIters.length; i++){
					SAMRecord sr_i = next(i);
					if(sr_i!=null) {
						ref_inds[i] = sr_i.getReferenceIndex(); //note this only advances if not returned;
						if(ref_inds[i]<minv){
							min_ind = i;
							minv = ref_inds[i];
						}
					}
				}
				current_sam_index = min_ind;
			}
			
			if(current_sam_index<0) return null;
			sr = this.currentVals[current_sam_index];
			if(sr!=null) this.currentIndex = sr.getReferenceIndex();
			cnts[current_sam_index]++;
			returned[current_sam_index] = true; 
			return sr;
			
		}
	}


//	private static final Logger LOG = LoggerFactory.getLogger(HTSErrorAnalysisCmd.class);

public static String getAnnotationsToInclude(String annotationType, boolean useExons){
	if(annotationType!=null) return annotationType;//.split(":");
	 if(!useExons) {
		 return  "gene:ncRNA:pseudogene:miRNA";	
	 }
	 else{
		 return "exon";
	 }
		
}
	public ViralTranscriptAnalysisCmd2() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("bamFile", null, "Name of bam file", true);
		addString("reference", null, "Name of reference genome", true);
		addString("annotation", null, "ORF annotation file or GFF file", false);
		addBoolean("useExons", true, "wehether to use exons");
		addString("readList", "", "List of reads", false);
		
			addString("annotType", null, "Type of annotation (only included if annotation is GFF file", false);
		addString("chroms", "all", "Restrict to these chroms, colon delimited", false);
		addString("resdir", "results"+System.currentTimeMillis(), "results directory");
		addString("GFF_features", "Name:description:ID:biotype:Parent", "GFF feature names");
		addBoolean("RNA", false, "If is direct RNA");
		addInt("maxReads", Integer.MAX_VALUE, "ORF annotation file");
		addInt("minClusterEntries",10,"threshold for consensus");
		addBoolean("tryComplementOnExtra", false, "look for negative strand matches on left over seqs");
		addBoolean("reAlignExtra", false, "whether to try realigning the extra sequence");
		addBoolean("combineOutput", false, "whether to combine output from different chroms");
		addString("pattern", null, "Pattern of read name, used for filtering");
		addString("span", "protein_coding", "Filtering span.  Use all not to filter.");
		addInt("qual", 0, "Minimum quality required");
		addInt("bin", 1, "Bin size for numerical hashing");
		addString("numBreaks", "none", "number of splices For  MSA");
		addInt("breakThresh", 1000, "Thresh for break points to match clusters.  If bigger than genome size then no break points");
		addBoolean("includeStart", true, "Whether to include start position in the cluster hash");

		addInt("coverageDepthThresh", 100, "Threshhold for writing base level depth information to h5 file");
		addInt("isoformDepthThresh", 10, "Threshhold for printing out all isoforms");
		addDouble("msaDepthThresh", 10, "Threshhold for running MSA per subcluster");
		addDouble("qualThresh", 20, "Quality thresh for leftover seqs");
		addString("doMSA", "false" , "Options: 5_3 or all or span=0 ");
		addString("msa_source", null , "indicates how to group msa, e.g 0,1,2;3,4,5 ");
		//addString("aligner", "kalign" , "Options: kalign3, poa, spoa, abpoa");
		addInt("startThresh", 100, "Threshold for having 5'");
		addInt("endThresh", 100, "Threshold for having 3'");
		addInt("maxThreads", 8, "max threads (only writing output uses threads at this stage)");
		addInt("extra_threshold", 200, "Threshold saving umatched 3'or 5'parts of reads");
		//addDouble("overlapThresh", 0.95, "Threshold for overlapping clusters");
	//	addBoolean("coexpression", false, "whether to calc coexperssion matrices (large memory for small bin size)");
		addBoolean("overwrite", false, "whether to overwrite existing results files");
		addBoolean("keepAlignment", false, "whether to keep alignment for MSA");
		addBoolean("attempt5rescue", true, "whether to attempt rescue of leader sequence if extra unmapped 5 read");
		addBoolean("attempt3rescue", true, "whether to attempt rescue of leader sequence if extra unmapped 5 read");
		addBoolean("writePolyA", false, "whether write reads with poly								qqq	A in middle");
		addBoolean("coronavirus", true, "whether to run in coronavirus mode (necessary to do breakpoint analysis, but takes more memory)");
		addStdHelp();
	}

	
 static void printParams(File resDir, String[] args1) throws IOException{
		File params_file = new File(resDir, "params.txt");
		PrintWriter pw = new PrintWriter(new FileWriter(params_file));
		Date d = new Date();
		
		 // System.getProperties().list(System.out);
		String mem = "Xmx="+(int)Math.ceil((double)Runtime.getRuntime().maxMemory()/(1e6))+"m";
		pw.println("User directory\n"+System.getProperty("user.dir"));
		pw.println("Date\n"+d.toString());
		pw.println("OS="+System.getProperty("os.name")+"\nJVM="+System.getProperty("java.version"));
		pw.println(mem);
		pw.println("Run commands");
		for(int i=0; i<args1.length; i+=1){
			pw.print(args1[i]+" ");
		}
		pw.close();
		
	}
	
 public static void run(CommandLine cmdLine, String[] bamFiles, String resDir,Map<String, JapsaAnnotation> anno,  Set<String> chrs) throws IOException{
		String reference = cmdLine.getStringVal("reference");
		int qual = cmdLine.getIntVal("qual");
		int bin = cmdLine.getIntVal("bin");
		int breakThresh = cmdLine.getIntVal("breakThresh");
		String pattern = cmdLine.getStringVal("pattern");
		String annotFile = cmdLine.getStringVal("annotation");
		String[] readList = cmdLine.getStringVal("readList").split(":");
		//String genesToInclude = cmdLine.getStringVal("genesToInclude");
		String  annotationType = getAnnotationsToInclude(cmdLine.getStringVal("annotType"), cmdLine.getBooleanVal("useExons"));
		boolean overwrite  = cmdLine.getBooleanVal("overwrite");
		int startThresh = cmdLine.getIntVal("startThresh");
		int endThresh = cmdLine.getIntVal("endThresh");
		int maxReads = cmdLine.getIntVal("maxReads");
		Annotation.enforceStrand = cmdLine.getBooleanVal("RNA");
		Outputs.executor=  cmdLine.getIntVal("maxThreads")==1 ? Executors.newSingleThreadExecutor():  Executors.newFixedThreadPool(cmdLine.getIntVal("maxThreads"));
//		Outputs.executor=  ;
		TranscriptUtils.qual_thresh = cmdLine.getDoubleVal("qualThresh");
	ViralTranscriptAnalysisCmd2.combineOutput = cmdLine.getBooleanVal("combineOutput");
		int isoformDepthThresh  = cmdLine.getIntVal("isoformDepthThresh");
		int coverageDepthThresh = cmdLine.getIntVal("coverageDepthThresh");
		IdentityProfile1.msaDepthThresh =(int) Math.floor(cmdLine.getDoubleVal("msaDepthThresh"));
	TranscriptUtils.extra_threshold = cmdLine.getIntVal("extra_threshold");
		TranscriptUtils.tryComplementOnExtra = cmdLine.getBooleanVal("tryComplementOnExtra");
		TranscriptUtils.reAlignExtra = cmdLine.getBooleanVal("reAlignExtra");
		TranscriptUtils.attempt5rescue = cmdLine.getBooleanVal("attempt5rescue");
		TranscriptUtils.attempt3rescue = cmdLine.getBooleanVal("attempt3rescue");
		Pattern patt = Pattern.compile(":");
		Outputs.numBreaks = patt.splitAsStream(cmdLine.getStringVal("numBreaks")=="none"  ? "" : cmdLine.getStringVal("numBreaks"))
		                            .map(Integer::valueOf)
		                            .collect(Collectors.toList());
		IdentityProfile1.includeStart = cmdLine.getBooleanVal("includeStart");
		GFFAnnotation.setGFFFeatureNames(cmdLine.getStringVal("GFF_features").split(":"));
		GFFAnnotation.span_only = cmdLine.getStringVal("span").equals("all") ?new ArrayList<String>() :   Arrays.asList(cmdLine.getStringVal("span").split(":"));
		boolean sorted = true;
		boolean coronavirus = cmdLine.getBooleanVal("coronavirus");
		String[] msaOpts = cmdLine.getStringVal("doMSA").split(":"); //e.g 5_3:sep or all:sep
		String msa_source = cmdLine.getStringVal("msa_source");
	//	String[] bamFiles =bamFile.split(":"); 
		int len =  bamFiles.length;
		if(msa_source!=null){
			String[] str = msa_source.split(";");
			outer: for(int j=0; j<bamFiles.length;j++){
				for(int i=0; i<str.length;i++){
						String[] str1 = str[i].split(",");
						for(int k=0; k<str1.length; k++){
							if(bamFiles[j].indexOf(str1[k])>=0){
								Outputs.msa_sources.put(j, i);
								continue outer;
							}
						}
					}
			}
			System.err.println("msa_sources "+Outputs.msa_sources);
		}else{
			
			for(int i=0; i<len; i++){
				Outputs.msa_sources.put(i, 0);
			}
		}
		for(int i=0; i<msaOpts.length; i++){
			if(msaOpts[i].startsWith("span=")) msaOpts[i] = msaOpts[i].split("=")[1];
		}
		Outputs.doMSA =Arrays.asList((msaOpts[0].equals("all") ?"5_3,no5_3,5_no3,no5_no3": msaOpts[0]).split(",")) ;
		Outputs.minClusterEntries = cmdLine.getIntVal("minClusterEntries");
		if(msaOpts[0].equals("false")){
			Outputs.doMSA = null;
		}
		Outputs.keepAlignment = cmdLine.getBooleanVal("keepAlignment");
		Outputs.keepinputFasta = Outputs.keepAlignment ;
		//Outputs.MSA_at_cluster = false;
		boolean calcBreaks=true;// = cmdLine.getBooleanVal("calcBreaks");// whether to calculate data for the break point heatmap, true for SARS_COV2
		boolean filterBy5_3 = false;// should be true for SARS_COV2
		boolean annotByBreakPosition = false;  // should be true for SARS_COV2
	//	String msa = cmdLine.getStringVal("aligner");
		//if(msa==null) throw new RuntimeException("!!");
	//	ErrorCorrection.msa = msa;
		Outputs.writePolyA = cmdLine.getBooleanVal("writePolyA");
		if(coronavirus){
		//	
			CigarCluster.recordDepthByPosition = true;
			System.err.println("running in coronavirus mode");
			calcBreaks  = true; 
			filterBy5_3 = true;
		//	Outputs.MSA_at_cluster = true;
			TranscriptUtils.checkAlign = true;
			TranscriptUtils.coronavirus = true;
			TranscriptUtils.extra_threshold1 =50;
			annotByBreakPosition = true;
			TranscriptUtils.writeAnnotP = true;
			CigarHash2.subclusterBasedOnStEnd = false;
		
		}else{
		//	TranscriptUtils.reAlignExtra = false;
		//	TranscriptUtils.findPolyA = false;
			CigarCluster.recordDepthByPosition = false;
			TranscriptUtils.coronavirus = false;
			TranscriptUtils.extra_threshold1 = 1000000;
			//Outputs.writeUnSplicedFastq = false;
			TranscriptUtils.checkAlign = false;
			TranscriptUtils.writeAnnotP = false;
			System.err.println("running in host mode");
			CigarHash2.subclusterBasedOnStEnd = false;
			calcBreaks = false;
		}
			errorAnalysis(bamFiles, reference, annotFile,readList,annotationType, 
				resDir,pattern, qual, bin, breakThresh, startThresh, endThresh,maxReads,  sorted , 
				calcBreaks, filterBy5_3, annotByBreakPosition, anno, chrs, overwrite, isoformDepthThresh, coverageDepthThresh);
	}
	public static void main(String[] args1) throws IOException, InterruptedException {
		long tme = System.currentTimeMillis();
		CommandLine cmdLine = new ViralTranscriptAnalysisCmd2();
		String[] args = cmdLine.stdParseLine(args1);
		String bamFile = cmdLine.getStringVal("bamFile");
		
		String resdir = cmdLine.getStringVal("resdir");
		String annot_file = cmdLine.getStringVal("annotation");
		String chroms= cmdLine.getStringVal("chroms");
		Set<String> chrs = null;
		if(!chroms.equals("all")){
			chrs = new HashSet<String>(Arrays.asList(chroms.split(":")));
		}
		boolean gff = annot_file !=null && annot_file.contains(".gff");
		Map<String, JapsaAnnotation> anno  = null;
		//boolean SARS = !gff;
		if(gff){
			String  annotationType = getAnnotationsToInclude(cmdLine.getStringVal("annotType"), cmdLine.getBooleanVal("useExons"));

			System.err.println("reading annotation");
			 anno =  GFFAnnotation.readAnno(annot_file, annotationType, chrs);
			// GFFAnnotation.read
			System.err.println("done reading annotation");
		}
		File resDir = new File(resdir);
		if(!resDir.exists())resDir.mkdir();
		printParams(resDir, args1);
		String[] bamFiles_ = bamFile.split(":");
		if(bamFile.equals("all") || bamFile.equals(".")){
			bamFiles_ = (new File("./")).list(new FilenameFilter(){

				@Override
				public boolean accept(File dir, String name) {
					return name.endsWith(".bam");
				}
				
			});
		}
		run(cmdLine, bamFiles_, resdir, anno,  chrs);
		Outputs.executor.shutdown();
		long time1 = System.currentTimeMillis();
		System.err.println((time1-tme)/(1000)+ " seconds");
		// paramEst(bamFile, reference, qual);
	}
public static boolean combineOutput = false;
	
	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String[] bamFiles_, String refFile, String annot_file, String[] readList,    String annotationType, String resdir, String pattern, int qual, int round, 
			int break_thresh, int startThresh, int endThresh, int max_reads,  boolean sorted,
			boolean calcBreaks , boolean filterBy5_3, boolean annotByBreakPosition,Map<String, JapsaAnnotation> anno, Set<String>chrToInclude, boolean overwrite,
			int writeIsoformDepthThresh, int writeCoverageDepthThresh) throws IOException {
		boolean cluster_reads = true;
		CigarHash2.round = round;
		
		IdentityProfile1.annotByBreakPosition = annotByBreakPosition;
		CigarHash.cluster_by_annotation =true;// cluster_by_annotation;
		TranscriptUtils.startThresh = startThresh;
		TranscriptUtils.endThresh = endThresh;
		IdentityProfile1.writeIsoformDepthThresh =writeIsoformDepthThresh;
		IdentityProfile1.writeCoverageDepthThresh =writeCoverageDepthThresh;
		File annotSummary = new File(resdir, "annotation.csv.gz");
		if(annotSummary.exists()) annotSummary.delete();
		PrintWriter annotation_pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(annotSummary, false))));
		boolean calcTree = false;
	//	String[] bamFiles_ = bamFiles.split(":");
	
		Collection<String>[] reads= null;
		if(readList.length>0 && readList[0].length()>0){
		 if( readList[0].indexOf("readToCluster.txt.gz")>=0){
			 Map<String, Collection<String>> map = new HashMap<String, Collection<String>>();
			 for(int i=0; i<readList.length; i++){
			 BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(new File(readList[i])))));
				String st = br.readLine();
				List<String> head = Arrays.asList(st.split("\\t"));
				int readind = head.indexOf("readID");
				int orfind = head.indexOf("ORFs");
				int chromind = head.indexOf("chrom");
				while((st = br.readLine())!=null){
					String[] str  = st.split("\\t");
					String readId = str[readind];
					String orfID = str[chromind]+";"+str[orfind];
					Collection<String> l= map.get(orfID) ;
					if(l==null) {
						map.put(orfID, new HashSet<String>());
						l= map.get(orfID) ;
					}
					l.add(readId);
				}
				br.close();
			 }
			readList = map.keySet().toArray(new String[0]);
			reads = new Collection[readList.length];
			for(int i=0; i<reads.length; i++){
				 reads[i]= map.get(readList[i]);
			}
		 }else{
			
		 reads= new Collection[readList.length];
		 for(int i=0; i<reads.length; i++){
			reads[i] = new HashSet<String>();
			BufferedReader br = new BufferedReader(new FileReader(new File(readList[i])));
			String st;
			while((st = br.readLine())!=null){
				String st_ = st.split("\\s+")[0];
				System.err.println(st_);
				reads[i].add(st_);
			}
			br.close();
		 }
		 }
		}
		
		TranscriptUtils.break_thresh = break_thresh;
		int len = bamFiles_.length;
		// len = 1;
		String[] in_nmes  = new String[len];
		ArrayList<Sequence> genomes = SequenceReader.readAll(refFile, Alphabet.DNA());

		// get the first chrom
		File resDir =  new File(resdir);
		if(!resDir.exists()) resDir.mkdir();
		else if(!resDir.isDirectory()) throw new RuntimeException(resDir+"should be a directory");
		///ArrayList<IdentityProfile1> profiles = new ArrayList<IdentityProfile1>();
		
		boolean gff = anno!=null;
		for (int ii = 0; ii < len; ii++) {
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
			in_nmes[ii] = bam.getName().split("\\.")[0];
		}
			
		
	//	genes_all_pw.close();
		IdentityProfile1 profile = null;
		Outputs outp = null;
		if(combineOutput)	outp = new Outputs(resDir,  in_nmes, overwrite, 0, true, CigarCluster.recordDepthByPosition); 
		final SAMRecordIterator[] samIters = new SAMRecordIterator[len];
		SamReader[] samReaders = new SamReader[len];
		outer1: for (int ii = 0; ii < len; ii++) {
			int source_index = ii;
			
			String bamFile = bamFiles_[ii];
			File bam = new File( bamFile);
		
			
			//https://www.programcreek.com/java-api-examples/?api=htsjdk.samtools.fastq.FastqWriter
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		//	SamReader samReader = null;// SamReaderFactory.makeDefault().open(new File(bamFile));

			if ("-".equals(bamFile))
				samReaders[ii] = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				samReaders[ii] = SamReaderFactory.makeDefault().open(bam);

			 samIters[ii] = samReaders[ii].iterator();
			// Read the reference genome
		}
		Set<Integer> chrom_indices_to_include = null;
		if(chrToInclude!=null){
			chrom_indices_to_include= new HashSet<Integer>();
			for(int i=0; i<genomes.size(); i++){
				if(chrToInclude.contains(genomes.get(i).getName())) chrom_indices_to_include.add(i);
			}
		}
		
		Iterator<SAMRecord> samIter= new CombinedIterator(samIters, max_reads,reads, chrom_indices_to_include);
			int currentIndex = 0;
			
			Sequence chr = genomes.get(currentIndex);
			int primelen = 500;//chr.length();
			Sequence chr5prime = TranscriptUtils.coronavirus ? chr.subSequence(0	, Math.min( primelen, chr.length())) : null;
			Sequence chr3prime = TranscriptUtils.coronavirus ? chr.subSequence(Math.max(0, chr.length()-primelen), chr.length()) : null;

			Set<String> doneChr = new HashSet<String>();
			
			long totReadBase = 0, totRefBase = 0;
			int numReads = 0;

			int numNotAligned = 0;
			int prev_src_index =-1;
			int numSecondary =0;
			outer: for (; samIter.hasNext() ; ) {
				SAMRecord sam= null;
				try{
					sam = samIter.next();
				}catch(Exception exc){
					exc.printStackTrace();
					
					
				}
				
				if(sam==null) break outer;
				int source_index = (Integer) sam.getAttribute(src_tag);
				int poolID = (Integer)  sam.getAttribute(pool_tag);
				if (pattern != null && (!sam.getReadName().contains(pattern)))
					continue;
				Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
				//String poolID = readSeq.get
				if (readSeq.length() <= 1) {
					//LOG.warn(sam.getReadName() +" ignored");
					
					// TODO: This might be secondary alignment, need to do something about it
					continue;
				}

				numReads++;
				if(sam.isSecondaryOrSupplementary()) {
					numSecondary++;
					continue;
				}
				if (sam.getReadUnmappedFlag()) {
					numNotAligned++;
					continue;
				}

				int flag = sam.getFlags();

				if (sam.getMappingQuality() < qual) {
					numNotAligned++;
					continue;
				}

				// int refPos = sam.getAlignmentStart() - 1;//convert to 0-based index
				int refIndex = sam.getReferenceIndex();
				
				// if move to another chrom, get that chrom
				if (refIndex != currentIndex) {
					if( profile!=null && currentIndex>=0){
						profile.printBreakPoints();
						profile.getConsensus();
						if(!combineOutput){
							outp.close();
							outp = null;
						}
						profile= null;
						
						doneChr.add(chr.getName());
						System.err.println("finished "+chr.getName());
						if(chrToInclude != null ){
							chrToInclude.remove(chr.getName());
							if(chrToInclude.size()==0) {
								System.err.println("finished sample "+in_nmes[source_index]);
								break outer;
							}
						}
					}
					currentIndex = refIndex;
					String prev_chrom = chr==null ? "null": chr.getName();
					if(currentIndex>=genomes.size() || currentIndex <0){
						System.err.println("for some reason the refIndex greater than genome size, so finishing here "+currentIndex+" "+genomes.size());
						break outer;
					}
					chr = genomes.get(currentIndex);
					chr5prime = TranscriptUtils.coronavirus ? chr.subSequence(0	, Math.min( primelen, chr.length())) : null;
					chr3prime = TranscriptUtils.coronavirus ? chr.subSequence(Math.max(0, chr.length()-primelen), chr.length()) : null;
					System.err.println("switch chrom "+prev_chrom+"  to "+chr.getName());
				}
				if(source_index!=prev_src_index){
					System.err.println("switch "+prev_src_index+"  to "+source_index+ " chr:"+chr.getName());
					 prev_src_index = source_index;
				}
				
				
					if(profile==null){
							if(doneChr.contains(chr.getName())){
								try{
									throw new RuntimeException("not sorted contains "+ chr.getName());
								}catch(Exception exc){
									exc.printStackTrace();
								}
							//	System.err.println("warning done"+chr.getName());
							}
							int seqlen = chr.length();
						//	if(!anno.containsKey(chr))
							Annotation annot = null;
						
							
							if(gff){
								JapsaAnnotation annot1 = anno.get(chr.getName());
								if(annot1==null) annot1 = anno.get("chr"+chr.getName());
								if(annot1==null) annot1 = anno.get(chr.getName().replaceAll("chr", ""));
								if(annot1==null){
									try{
										annot = new EmptyAnnotation(chr.getName(), chr.getDesc(), seqlen, annotation_pw);
									System.err.println("no annotation for  "+ chr.getName());
									}catch(Exception exc){
										exc.printStackTrace();
									}
								}else{
									
									annot = new GFFAnnotation(annot1, seqlen, annotation_pw, len);
									
								}
							}else{
								annot = annot_file == null ? new EmptyAnnotation(chr.getName(), chr.getDesc(), seqlen, annotation_pw) : 
									new Annotation(new File(annot_file), currentIndex+"", seqlen, len);
							}
						//	pw.close();
						if(!combineOutput)	outp = new Outputs(resDir,  in_nmes, overwrite, currentIndex, true, CigarCluster.recordDepthByPosition); 
							boolean calcBreaks1 = calcBreaks && break_thresh < seqlen;
							profile = new IdentityProfile1(chr, outp,  in_nmes, startThresh, endThresh, annot, calcBreaks1, chr.getName(), currentIndex);
					}
				
					try{
						String pool = readList==null || poolID<0 ? "" : (readList[poolID]+"|");
						TranscriptUtils.identity1(chr, chr5prime,chr3prime, readSeq, sam, profile, source_index, cluster_reads, chr.length(), pool);
					}catch(NumberFormatException exc){
						System.err.println(readSeq.getName());
						exc.printStackTrace();
					}
			}//samIter.hasNext()
			for(int ii=0; ii<samReaders.length; ii++){
			samReaders[ii].close();
			}
			
		
		
			if(profile!=null){
				profile.printBreakPoints();
				profile.getConsensus();
				
			}
		
		if(outp!=null) outp.close();
	
		if(annotation_pw!=null) annotation_pw.close();
		
	
	}
	
}
