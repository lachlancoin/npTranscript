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
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.Annotation;
import npTranscript.cluster.CigarHash;
import npTranscript.cluster.CigarHash2;
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
private static final class CombinedIterator implements Iterator<SAMRecord> {
		private final SAMRecordIterator[] samIters;
		int current_sam_index =0;
		int currentIndex =0; //relates to chromosome
		private final SAMRecord[] currentVals;
		private final boolean[] returned;
		private final int[] cnts ;
		int max;
		Set<String> readList ;
		Set<Integer> chrs;
		private CombinedIterator(SAMRecordIterator[] samIters, int max, Set<String>readList, Set<Integer> chrs) {
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
						if(readList!=null && readList.contains(sr.getReadName())) break inner;
						if(readList==null &&  chrs!=null && chrs.contains(sr.getReferenceIndex())) break inner;
						if(readList==null  && chrs==null) break inner;
					}
					if(sr!=null) {
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

	
	public ViralTranscriptAnalysisCmd2() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("bamFile", null, "Name of bam file", true);
		addString("reference", null, "Name of reference genome", true);
		addString("annotation", null, "ORF annotation file or GFF file", true);
		addString("readList", null, "List of reads", false);
		String all_types = "gene:ncRNA_gene:pseudogene";			
			addString("type", all_types, "Type of annotation (only included if annotation is GFF file", false);
		addString("chroms", "all", "Restrict to these chroms, colon delimited", false);
		addString("resdir", "results"+System.currentTimeMillis(), "results directory");
		addInt("maxReads", Integer.MAX_VALUE, "ORF annotation file");
		addString("pattern", null, "Pattern of read name, used for filtering");
		addInt("qual", 0, "Minimum quality required");
		addInt("bin", 1, "Bin size for numerical hashing");
		addInt("breakThresh", 10, "Thresh for break points to match clusters");
		addInt("coverageDepthThresh", 100, "Threshhold for writing base level depth information to h5 file");
		addInt("isoformDepthThresh", 10, "Threshhold for printing out all isoforms");
		addDouble("msaDepthThresh", 10, "Threshhold for running MSA per subcluster");
		addString("doMSA", "false" , "Options: combined, separate or false");
		
		addInt("startThresh", 100, "Threshold for having 5'");
		addInt("endThresh", 100, "Threshold for having 3'");
		addInt("extra_threshold", 200, "Threshold saving umatched 3'or 5'parts of reads");
		//addDouble("overlapThresh", 0.95, "Threshold for overlapping clusters");
	//	addBoolean("coexpression", false, "whether to calc coexperssion matrices (large memory for small bin size)");
		addBoolean("overwrite", false, "whether to overwrite existing results files");
		addBoolean("keepAlignment", false, "whether to keep alignment for MSA");
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
	
 public static void run(CommandLine cmdLine, String bamFile, String resDir,Map<String, JapsaAnnotation> anno,  Set<String> chrs) throws IOException{
		String reference = cmdLine.getStringVal("reference");
		int qual = cmdLine.getIntVal("qual");
		int bin = cmdLine.getIntVal("bin");
		int breakThresh = cmdLine.getIntVal("breakThresh");
		String pattern = cmdLine.getStringVal("pattern");
		String annotFile = cmdLine.getStringVal("annotation");
		String readList = cmdLine.getStringVal("readList");
		//String genesToInclude = cmdLine.getStringVal("genesToInclude");
		String  annotationType = cmdLine.getStringVal("type");
		boolean overwrite  = cmdLine.getBooleanVal("overwrite");
		int startThresh = cmdLine.getIntVal("startThresh");
		int endThresh = cmdLine.getIntVal("endThresh");
		int maxReads = cmdLine.getIntVal("maxReads");
	
		int isoformDepthThresh  = cmdLine.getIntVal("isoformDepthThresh");
		int coverageDepthThresh = cmdLine.getIntVal("coverageDepthThresh");
		IdentityProfile1.msaDepthThresh =(int) Math.floor(cmdLine.getDoubleVal("msaDepthThresh"));
	TranscriptUtils.extra_threshold = cmdLine.getIntVal("extra_threshold");
		
		boolean sorted = true;
		boolean coronavirus = cmdLine.getBooleanVal("coronavirus");
		String[] msaOpts = cmdLine.getStringVal("doMSA").split(":"); //e.g 5_3:sep or all:sep
		
		Outputs.doMSA =Arrays.asList((msaOpts[0].equals("all") ?"5_3:no5_3:5_no3:no5_no3": msaOpts[0]).split(":")) ;
		if(msaOpts[0].equals("false")){
			Outputs.doMSA = null;
		}else{
			Outputs.mergeSourceClusters = msaOpts.length==1 || !msaOpts[1].startsWith("sep");
		}
		Outputs.keepAlignment = cmdLine.getBooleanVal("keepAlignment");
		Outputs.keepinputFasta = Outputs.keepAlignment ;
		//Outputs.MSA_at_cluster = false;
		boolean calcBreaks = false;// whether to calculate data for the break point heatmap, true for SARS_COV2
		boolean filterBy5_3 = false;// should be true for SARS_COV2
		boolean annotByBreakPosition = false;  // should be true for SARS_COV2
	
		
		if(coronavirus){
			System.err.println("running in coronavirus mode");
			calcBreaks  = true; 
			filterBy5_3 = true;
		//	Outputs.MSA_at_cluster = true;
			TranscriptUtils.checkAlign = true;
			annotByBreakPosition = true;
			CigarHash2.subclusterBasedOnStEnd = true;
		}else{
			TranscriptUtils.checkAlign = false;
			System.err.println("running in host mode");
			CigarHash2.subclusterBasedOnStEnd = false;
		}
			errorAnalysis(bamFile, reference, annotFile,readList,annotationType, 
				resDir,pattern, qual, bin, breakThresh, startThresh, endThresh,maxReads,  sorted , 
				calcBreaks, filterBy5_3, annotByBreakPosition, anno, chrs, overwrite, isoformDepthThresh, coverageDepthThresh);
	}
	public static void main(String[] args1) throws IOException, InterruptedException {
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
		boolean gff = annot_file.contains(".gff");
		Map<String, JapsaAnnotation> anno  = null;
		//boolean SARS = !gff;
		if(gff){
			String  annotationType = cmdLine.getStringVal("type");
			System.err.println("reading annotation");
			 anno =  GFFAnnotation.readAnno(annot_file, annotationType, chrs);
			// GFFAnnotation.read
			System.err.println("done reading annotation");
		}
		File resDir = new File(resdir);
		if(!resDir.exists())resDir.mkdir();
		printParams(resDir, args1);
		run(cmdLine, bamFile, resdir, anno,  chrs);
		
		
		// paramEst(bamFile, reference, qual);
	}

	
	/**
	 * Error analysis of a bam file. Assume it has been sorted
	 */
	static void errorAnalysis(String bamFiles, String refFile, String annot_file, String readList,    String annotationType, String resdir, String pattern, int qual, int round, 
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

		boolean calcTree = false;
		String[] bamFiles_ = bamFiles.split(":");
		Set<String> reads = null;
		if(readList!=null){
			reads = new HashSet<String>();
			BufferedReader br = new BufferedReader(new FileReader(new File(readList)));
			String st;
			while((st = br.readLine())!=null){
				String st_ = st.split("\\s+")[0];
				System.err.println(st_);
				reads.add(st_);
			}
			br.close();
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
			int primelen = 1000;
			Sequence chr3prime = chr.subSequence(chr.length()-primelen, chr.length());
			Sequence chr5prime = chr.subSequence(0	, primelen);

			Set<String> doneChr = new HashSet<String>();
			
			long totReadBase = 0, totRefBase = 0;
			int numReads = 0;

			int numNotAligned = 0;
			int prev_src_index =-1;
			
			outer: for (; samIter.hasNext() ; ) {
				SAMRecord sam= null;
				try{
					sam = samIter.next();
				}catch(Exception exc){
					exc.printStackTrace();
					
					
				}
				if(sam==null) break outer;
				int source_index = (Integer) sam.getAttribute(src_tag);
			//	sam.getBaseQualityString();
				//System.err.println(source_index);
				if (pattern != null && (!sam.getReadName().contains(pattern)))
					continue;

				// make the read seq
				Sequence readSeq = new Sequence(Alphabet.DNA(), sam.getReadString(), sam.getReadName());
				if (readSeq.length() <= 1) {
					// LOG.warn(sam.getReadName() +" ignored");
					// TODO: This might be secondary alignment, need to do something about it
					continue;
				}

				numReads++;

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
						outp.close();
						profile= null;
						outp = null;
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
					chr = genomes.get(currentIndex);
					 chr3prime = chr.subSequence(chr.length()-primelen, chr.length());
					chr5prime = chr.subSequence(0	, primelen);
				//	outp.updateChromIndex(currentIndex);
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
								if(annot1==null){
									try{
										throw new RuntimeException("no annotation for  "+ chr.getName());
									}catch(Exception exc){
										exc.printStackTrace();
									}
								}else{
									annot = new GFFAnnotation(annot1, seqlen);
								}
							}else{
								annot = new Annotation(new File(annot_file), currentIndex+"", seqlen);
							}
							outp = new Outputs(resDir,  in_nmes, overwrite, currentIndex, true, true, seqlen); 
							profile = new IdentityProfile1(chr, outp,  in_nmes, startThresh, endThresh, annot, calcBreaks, chr.getName(), currentIndex);
					}
					try{
						TranscriptUtils.identity1(chr, chr5prime, chr3prime, readSeq, sam, profile, source_index, cluster_reads, chr.length());
					}catch(NumberFormatException exc){
						System.err.println(readSeq.getName());
						exc.printStackTrace();
					}
			}
			for(int ii=0; ii<samReaders.length; ii++){
			samReaders[ii].close();
			}
			
		
		
			if(profile!=null){
				profile.printBreakPoints();
				profile.getConsensus();
				
			}
		
		if(outp!=null) outp.close();
		File[] reads_files = resDir.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				// TODO Auto-generated method stub
				return pathname.getName().contains("readToCluster.txt");
			}
			
		});
		
		for(int i=0; i<reads_files.length; i++){
			System.err.println(reads_files[i]);
			File read_fle = reads_files[i];
			final String prefix = read_fle.getName().split("\\.")[0]+".";
			System.err.println(prefix);
			File[] transcript_files = resDir.listFiles(new FileFilter(){

				@Override
				public boolean accept(File pathname) {
					// TODO Auto-generated method stub
					return pathname.getName().contains(".transcripts.txt") && pathname.getName().startsWith(prefix);
				}
				
			});
			/*if(transcript_files.length==1){
				ExtractClusterCmd.cluster(resDir,transcript_files[0],reads_files[i], filterBy5_3);
			}*/
		}
	
	}

}

/*
 * RST* ----------------------------------------------------------
 * jsa.hts.errorAnalysis*: Error analysis of sequencing data
 * ----------------------------------------------------------
 * 
 * jsa.hts.errorAnalysis* assesses the error profile of sequencing data by
 * getting the numbers of errors (mismatches, indels etc) from a bam file.
 * Obviously, it does not distinguish sequencing errors from mutations, and
 * hence consider mutations as errors. It is best to use with the bam file from
 * aligning sequencing reads to a reliable assembly of the sample.
 * 
 * <usage>
 * 
 * RST
 */