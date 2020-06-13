/*****************************************************************************
 * Copyright (c)Lachlan Coin University of Melbourne, All rights reserved.   *
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
 * 05/05/2020 - Chenxi Zhou                                       
 ****************************************************************************/

package npTranscript.run;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.SequenceUtil;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.util.ConcurrencyTools;

/**
 * @author Chenxi Zhou
 *
 */

@Deployable(scriptName = "npChimera.run", scriptDesc = "Detection of chimeric reads")
public class ViralChimericReadsAnalysisCmd extends CommandLine {

	private static boolean IN_DEBUG_MODE;
	
	public ViralChimericReadsAnalysisCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addBoolean("getClippedReads", false, "Retrieve reads unmapped or mapped with clippings");
		addBoolean("getLeftover", false, "Retrieve leftover sequences");
		addBoolean("postAnalysis", false, "Post analysis for chimeric reads");
		addString("bamFile", null, "Name of BAM file", true);
		addString("outPrefix", null, "Prefix for output files", true);
		addString("reference", null, "Reference file in FASTA format", true);
		addInt("minReadLen", 500, "Minimum read length");
		addInt("hitThresh", 200, "Length threshold for hits on the genome");
		addInt("clipThresh", 20, "Length threshold for clipping");
		addInt("leftOverThresh", 200, "Length threshold for left over sequences");
		addInt("qualThresh", 0, "Quality score threshold for clipping");
		addInt("flankSize", 20, "Flanking size for the reference sequence at fusion point");
		addInt("minClippedLen", 200, "Minimum length of clippings (only for 'getClippedReads')");
		addBoolean("debug", false, "Running in debug mode.");
		addStdHelp();
	}
	
	public static void main(String[] args) {
		CommandLine cmdLine = new ViralChimericReadsAnalysisCmd();
		cmdLine.stdParseLine(args);
		getClippedReads = cmdLine.getBooleanVal("getClippedReads");
		getLeftover = cmdLine.getBooleanVal("getLeftover");
		postAnalysis = cmdLine.getBooleanVal("postAnalysis");
		int comds = 0;
		comds += getClippedReads ? 1 : 0;
		comds += getLeftover ? 1 : 0;
		comds += postAnalysis ? 1 : 0;
		if(comds == 0)
			throw new RuntimeException("Please specify one option: --getClippedReads, --getLeftover, --postAnalysis"); 
		if(comds > 1)
			throw new RuntimeException("Options --getClippedReads, --getLeftover and --postAnalysis are mutually exclusive"); 
		
		String bamFile = cmdLine.getStringVal("bamFile");
		String outPrefix = cmdLine.getStringVal("outPrefix");
		String reference = cmdLine.getStringVal("reference");
		min_readLen = cmdLine.getIntVal("minReadLen");
		hit_thresh = cmdLine.getIntVal("hitThresh");
		clip_thresh = cmdLine.getIntVal("clipThresh");
		leftover_thresh = cmdLine.getIntVal("leftOverThresh");
		qual_thresh = cmdLine.getIntVal("qualThresh");
		IN_DEBUG_MODE = cmdLine.getBooleanVal("debug");
		min_clipped = cmdLine.getIntVal("minClippedLen");
		
		if(clip_thresh>leftover_thresh) {
			System.err.println("WRANING: leftOverThresh is smaller than clipThresh. "
					+ "Set leftOverThresh="+clip_thresh);
			leftover_thresh = clip_thresh;
		}
		flank_size = cmdLine.getIntVal("flankSize");
		run(bamFile.split(":"), reference, outPrefix);
	}
	
	private static boolean getClippedReads;
	private static boolean getLeftover;
	private static boolean postAnalysis;
	private static int min_readLen;
	private static int hit_thresh;
	private static int leftover_thresh;
	private static int clip_thresh;
	private static double qual_thresh;
	private static int flank_size;
	private static int min_clipped;
	
	private static class Block {
		private String str_id;
		private int s1;
		private int s2;
		private int q1;
		private int q2;
		private String strand;

		public Block(String str_id, int s1, int s2, int q1, int q2, String strand) {
			this.str_id = str_id;
			this.s1 = s1;
			this.s2 = s2;
			this.q1 = q1;
			this.q2 = q2;
			this.strand = strand;
		}
	}
	
	public static void run(String[] files, String reference, String out_prefix) {
		// TODO Auto-generated method stub
		String out_clips = getClippedReads ? out_prefix+".fastq" : null;
		String out_fq = getLeftover ? out_prefix+"_complete.fastq" : null;
		String out_leftover = getLeftover ? out_prefix+"_leftover.fastq" : null;
		String out_bam = getLeftover ? out_prefix+".bam" : null;
		String out_qual = (IN_DEBUG_MODE && getLeftover) ? out_prefix+"_qual.txt" : null;
		String out_aln = postAnalysis ? out_prefix+".aln" : null;
		if((out_fq != null && (new File(out_fq)).exists()) || 
				out_leftover != null && (new File(out_leftover)).exists() ||
				out_bam != null && (new File(out_bam)).exists() ||
				out_clips != null && (new File(out_clips)).exists() || 
				out_qual != null && (new File(out_qual)).exists() || 
				out_aln != null && (new File(out_aln)).exists())
			throw new RuntimeException("Output file exsits!!!");
		
		try {
			BufferedWriter bw_complete = out_fq == null ? null : new BufferedWriter(new FileWriter(new File(out_fq)));
			BufferedWriter bw_leftover = out_leftover == null ? null : new BufferedWriter(new FileWriter(new File(out_leftover)));
			BufferedWriter bw_clips = out_clips == null ? null : new BufferedWriter(new FileWriter(new File(out_clips)));
			BufferedWriter bw_qual = out_qual == null ? null : new BufferedWriter(new FileWriter(new File(out_qual)));
			BufferedWriter bw_aln = out_aln == null ? null : new BufferedWriter(new FileWriter(new File(out_aln)));
			IndexedFastaSequenceFile dict = getClippedReads ? null : new IndexedFastaSequenceFile(new File(reference));
			SAMFileWriter samWriter = null;
			for(String f : files) {
				List<String> bam_files = new ArrayList<>();
				File file = new File(f);
				if(!file.exists()) {
					if(bw_complete != null) bw_complete.close(); 
					if(bw_leftover != null) bw_leftover.close(); 
					if(bw_clips != null) bw_clips.close(); 
					if(bw_qual != null) bw_qual.close(); 
					if(bw_aln != null) bw_aln.close(); 
					if(dict != null) dict.close(); 
					if(samWriter != null) samWriter.close(); 
					throw new RuntimeException("!!!");
				}
				
				if(file.isFile()) {
					bam_files.add(f);
				} else if(file.isDirectory()) {
					File dir = new File(f);
					File[] fs = dir.listFiles(new FilenameFilter() {
						public boolean accept(File dir, String name) {
							return name.toLowerCase().endsWith(".bam");
						}
					});
					for(File s : fs)
						bam_files.add(s.getAbsolutePath());
				}
				
				for(String bam_file : bam_files) {
					SamReader samReader = SamReaderFactory.make()
							.validationStringency(ValidationStringency.LENIENT)
							.samRecordFactory(DefaultSAMRecordFactory.getInstance())
							.open(new File(bam_file));
					if(out_bam != null && samWriter == null)
						samWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReader.getFileHeader(),
								true, new File(out_bam));
					SAMRecordIterator iterator = samReader.iterator();
					if(qual_thresh == 0.0D) {
						qual_thresh = getLeftover ? calcQSMeanAndPercentiles(bam_file)[26] : 0.0D;
						System.err.println("#quality score q25: " + qual_thresh);
					}
					
					List<SAMRecord> records = new ArrayList<>();
					Map<String, List<Block>> blocks = new HashMap<>();
					Map<String, RangeSet<Integer>> ranges = new HashMap<>();
					SAMRecord record = iterator.hasNext() ? (SAMRecord)iterator.next() : null;
					String read_id = (record == null) ? null : record.getReadName();
					while(record != null) {
						if(record.getReadUnmappedFlag()) {
							if(getClippedReads && record.getReadLength() >= min_readLen)
								bw_clips.write("@" + read_id + "\n" + record.getReadString() + "\n+\n" + record.getBaseQualityString() + "\n"); 
							record = iterator.hasNext() ? (SAMRecord)iterator.next() : null;
							read_id = (record == null) ? null : record.getReadName();
							continue;
						} 
						records.clear();
						records.add(record);
						while(true) {
							if((record = iterator.hasNext() ? iterator.next() : null) == null || 
									!record.getReadName().equals(read_id))
								break; 
							records.add(record);
						} 
						blocks.clear();
						ranges.clear();
						int readLen = 0;
						String sequence = null;
						String qual = null;
						byte[] phredQs = null;
						String tagLO = null;
						String tagAP = null;
						String tagFR = null;
						String tagFQ = null;
						String tagRS = null;
						for(SAMRecord r : records) {
							if(!r.isSecondaryOrSupplementary()) {
								readLen = r.getReadLength();
								sequence = r.getReadString();
								qual = r.getBaseQualityString();
								phredQs = r.getBaseQualities();
								if(r.getReadNegativeStrandFlag()) {
									sequence = SequenceUtil.reverseComplement(sequence);
									qual = (new StringBuilder(qual)).reverse().toString();
									reverseArray(phredQs);
								} 
								tagLO = r.getStringAttribute("LO");
								tagAP = r.getStringAttribute("AP");
								tagFR = r.getStringAttribute("FR");
								tagFQ = r.getStringAttribute("FQ");
								tagRS = r.getStringAttribute("RS");
							} 
							int hclip = 0;
							CigarElement c = r.getCigar().getFirstCigarElement();
							if(c.getOperator() == CigarOperator.HARD_CLIP)
								hclip = c.getLength(); 
							int s1 = r.getAlignmentStart();
							int s2 = r.getAlignmentEnd();
							int q1 = r.getReadPositionAtReferencePosition(s1) + hclip;
							int q2 = r.getReadPositionAtReferencePosition(s2) + hclip;
							if(r.getReadNegativeStrandFlag()) {
								int tmp = q1;
								q1 = readLen - q2 + 1;
								q2 = readLen - tmp + 1;
							} 
							String str = r.getReferenceName();
							ranges.putIfAbsent(str, TreeRangeSet.create());
							blocks.putIfAbsent(str, new ArrayList<>());
							ranges.get(str).add(Range.closed(q1, q2).canonical(DiscreteDomain.integers()));
							blocks.get(str).add(new Block(r.getReferenceName(), s1, s2, q1, q2, r.getReadNegativeStrandFlag() ? "r" : "f"));
						} 
						if(readLen < min_readLen && !postAnalysis) {
							if(record != null)
								read_id = record.getReadName(); 
							continue;
						} 
						String refId = null;
						int alnLen = 0;
						for(String ref : ranges.keySet()) {
							int len = 0;
							for(Range<Integer> r : ranges.get(ref).asRanges())
								len += r.upperEndpoint() - r.lowerEndpoint(); 
							if(len > alnLen) {
								alnLen = len;
								refId = ref;
							} 
						} 
						if(alnLen < hit_thresh) {
							if(getClippedReads && readLen - alnLen >= min_clipped)
								bw_clips.write("@" + read_id + "\n" + sequence + "\n+\n" + qual + "\n"); 
							if(record != null)
								read_id = record.getReadName(); 
							continue;
						} 
						List<Block> block = blocks.get(refId);
						Collections.sort(block, new Comparator<Block>() {
							public int compare(ViralChimericReadsAnalysisCmd.Block o1, ViralChimericReadsAnalysisCmd.Block o2) {
								return Integer.compare(o1.q1, o2.q1);
							}
						});
						Block blk = block.get(0);
						boolean isPositiveStrand = blk.strand.equals("f");
						int qstart = blk.q1;
						int qend = blk.q2;
						int sstart = isPositiveStrand ? blk.s1 : blk.s2;
						int send = isPositiveStrand ? blk.s2 : blk.s1;
						boolean badAlignment = false;
						for(int j = 1; j < block.size(); j++) {
							blk = block.get(j);
							int olap = isPositiveStrand ? (send - blk.s1) : (blk.s2 - send);
							if(!blk.str_id.equals(refId) || 
									isPositiveStrand != blk.strand.equals("f") || 
									Math.abs(qend - blk.q1) > 20 || 
									olap > 20) {
								badAlignment = true;
								break;
							} 
							qend = blk.q2;
							send = isPositiveStrand ? blk.s2 : blk.s1;
						} 
						
						// 1. check if the alignment chunks are consistent (for supplementary alignments)
						//     - alignment strand need to be consistent for all alignment chunks
						//     - no gaps or overlap over 20bp on read sequence between multiple alignment chunks
						//     - no overlap over 20bp on reference between multiple alignment chunks (gaps are allowed)
						// 2. fusion in only one side
						//     - if the clip at 5'- or 3'-end is over 20bp, will check if it's polyA tail or not
						//       if not will consider as a clip
						//     - clip is only allowed in one side
						// 3. output of leftover sequences
						//     - if the clip sequence is over 200bp, will output as leftover sequences
						//       leftover sequence will be only at 5'-end or 3'-end (otherwise will be filtered out by #2)
						//     - the reference sequence flanking 20bp of the upstream and downstream at the fusion point
						//       will be always 41bp by default - padded with 'N' if flanking sequence is too short
						// 4. the FASTQ header of the leftover sequences will contain four information fields
						//     - add '-y' option for alignment with minimap2
						//     - SAMRecord getStringAttribute() could be used to retrieve these information
                        //       record.getStringAttribute("LO") for the position of the leftover sequence
						//       record.getStringAttribute("AP") for the alignment position of the original sequence
						//       record.getStringAttribute("FR") for the flanking sequence on the reference sequence
						//       record.getStringAttribute("FQ") for the flanking sequence on the query/read sequence
						//       record.getStringAttribute("RS") for the reference sequence name
						
						boolean leftOverAt5End = (qstart > clip_thresh && !isPolyATail(sequence.substring(0, qstart - 1)));
						boolean leftOverAt3End = (qend <= readLen - clip_thresh && !isPolyATail(sequence.substring(qend)));
						if(getClippedReads) {
							if(badAlignment || leftOverAt5End || leftOverAt3End)
								bw_clips.write("@" + read_id + "\n" + sequence + "\n+\n" + qual + "\n"); 
							if(record != null)
								read_id = record.getReadName(); 
							continue;
						} 
						if(getLeftover) {
							if(badAlignment || leftOverAt5End == leftOverAt3End) {
								if(record != null)
									read_id = record.getReadName(); 
								continue;
							} 
							
							// build describe line
							StringBuilder descBuilder = new StringBuilder("AP:Z:");
							for(Block block1 : block) {
								descBuilder.append(block1.q1);
								descBuilder.append(",");
								descBuilder.append(block1.q2);
								descBuilder.append(",");
								descBuilder.append(block1.s1);
								descBuilder.append(",");
								descBuilder.append(block1.s2);
								descBuilder.append(",");
								descBuilder.append(block1.strand);
								descBuilder.append(";");
							} 
							descBuilder.setLength(descBuilder.length() - 1);
							String desc = descBuilder.toString();
							String refSequence = dict.getSequence(refId).getBaseString();
							boolean isChimera = false;
							if(qstart > leftover_thresh && leftOverAt5End) {
								String seqstr = sequence.substring(0, qstart - 1);
								String qualstr = qual.substring(0, qstart - 1);
								// check quality score in downstream region
								if(median(phredQs, Math.max(0, qstart - 31), Math.min(30, qstart - 1)) >= qual_thresh && 
										median(phredQs, Math.max(0, qstart - 51), Math.min(50, qstart - 1)) >= qual_thresh && 
										median(phredQs, Math.max(0, qstart - 71), Math.min(70, qstart - 1)) >= qual_thresh && 
										median(phredQs, Math.max(0, qstart - 101), Math.min(100, qstart - 1)) >= qual_thresh && 
										median(phredQs, 0, qstart - 1) >= qual_thresh) {
									String fqstr = getFlankSequence(sequence, qstart - 1, flank_size);
									String frstr = getFlankSequence(refSequence, sstart - 1, flank_size);
									if(!isPositiveStrand)
										frstr = SequenceUtil.reverseComplement(frstr); 
									bw_leftover.write("@" + read_id + "\tLO:Z:" + 1 + "," + (qstart - 1) + "," + readLen + "\t" + desc + "\t" + 
											"FR:Z:" + frstr + "," + sstart + "," + (isPositiveStrand ? "f" : "r") + "\t" + 
											"FQ:Z:" + fqstr + "," + qstart + ",3'\t" + 
											"RS:Z:" + refId + "\n"); // match at 3'-end
									bw_leftover.write(seqstr + "\n");
									bw_leftover.write("+\n");
									bw_leftover.write(qualstr + "\n");
									if(IN_DEBUG_MODE) {
										String fqQualStr = getFlankSequence(qual, qstart - 1, 300, '!');
										//bw_qual.write(read_id+"|"+fqstr+"|"+fqQualStr+"|"+qstart+"|3'\t")
										bw_qual.write(read_id);
										byte[] qs = SAMUtils.fastqToPhred(fqQualStr);
										//reverse quality str to make the flanking position at the 3'-end
										reverseArray(qs);
										for(byte b : qs ) bw_qual.write("\t" + b);
										bw_qual.write("\n");
									} 
									isChimera = true;
								} 
							} 
							if(qend <= readLen - leftover_thresh && leftOverAt3End) {
								String seqstr = sequence.substring(qend);
								String qualstr = qual.substring(qend);
								int leftLen = readLen - qend;
								// check quality score in downstream region
								if(median(phredQs, qend, Math.min(30, leftLen)) >= qual_thresh && 
										median(phredQs, qend, Math.min(50, leftLen)) >= qual_thresh && 
										median(phredQs, qend, Math.min(70, leftLen)) >= qual_thresh && 
										median(phredQs, qend, Math.min(100, leftLen)) >= qual_thresh && 
										median(phredQs, qend, leftLen) >= qual_thresh) {
									String fqstr = getFlankSequence(sequence, qend - 1, flank_size);
									String frstr = getFlankSequence(refSequence, send - 1, flank_size);
									if(!isPositiveStrand)
										frstr = SequenceUtil.reverseComplement(frstr); 
									bw_leftover.write("@" + read_id + "\tLO:Z:" + (qend + 1) + "," + readLen + "," + readLen + "\t" + desc + "\t" + 
											"FR:Z:" + frstr + "," + send + "," + (isPositiveStrand ? "f" : "r") + "\t" + 
											"FQ:Z:" + fqstr + "," + qend + ",5'\t" + 
											"RS:Z:" + refId + "\n"); // match at 5'-end
									bw_leftover.write(seqstr + "\n");
									bw_leftover.write("+\n");
									bw_leftover.write(qualstr + "\n");
									if(IN_DEBUG_MODE) {
										String fqQualStr = getFlankSequence(qual, qend - 1, 300, '!');
										//bw_qual.write(read_id+"|"+fqstr+"|"+fqQualStr+"|"+qend+"|5'");
										bw_qual.write(read_id);
										byte[] qs = SAMUtils.fastqToPhred(fqQualStr);
										for(byte b : qs) bw_qual.write("\t" + b);
										bw_qual.write("\n");
									} 
									isChimera = true;
								} 
							} 
							if(isChimera) {
								bw_complete.write("@" + read_id + "\t" + desc + "\t" + "RS:Z:" + refId + "\n");
								bw_complete.write(sequence + "\n");
								bw_complete.write("+\n");
								bw_complete.write(qual + "\n");
								for(SAMRecord r : records)
									samWriter.addAlignment(r); 
							} 
						} else if(postAnalysis) {
							if(badAlignment) {
								if(record != null)
									read_id = record.getReadName(); 
								continue;
							} 
							String[] s = tagLO.split(",");
							int rlen = Integer.parseInt(s[2]);
							s = tagAP.split(",");
							int qs = Integer.parseInt(s[0]);
							int qe = Integer.parseInt(s[1]);
							s = tagFR.split(",");
							String refSeq = s[0];
							s = tagFQ.split(",");
							String qrySeq = s[0];
							boolean seqAt5End = !s[2].equals("5'");
							String refSequence = dict.getSequence(refId).getBaseString();
							String[] seqs = new String[3];
							if(seqAt5End) {
								String fqstr = getFlankSequence(sequence, qend - 1, flank_size);
								String frstr = getFlankSequence(refSequence, send - 1, flank_size);
								if(!isPositiveStrand)
									frstr = SequenceUtil.reverseComplement(frstr); 
								String gapSeq = sequence.substring(qend);
								String qstr = fqstr.substring(0, fqstr.length() / 2 + 1) + 
										gapSeq + 
										qrySeq.substring(qrySeq.length() / 2);
								bw_aln.write("@" + read_id + " AP:Z:" + (qend - flank_size) + "," + qend + "," + qs + "," + (qs + qrySeq.length() / 2) + " RL:i:" + rlen + "\n");
								bw_aln.write("#AP:Z:" + qstart + "," + qend + "," + sstart + "," + send + "," + (isPositiveStrand ? "f" : "r") + 
										" FR:Z:" + frstr + "," + send + "," + (isPositiveStrand ? "f" : "r") + 
										" FQ:Z:" + fqstr + "," + qend + ",5'" + 
										" RS:Z:" + refId + "\n");
								bw_aln.write("#AP:Z:" + tagAP + " FR:Z:" + tagFR + " FQ:Z:" + tagFQ + " RS:Z:" + tagRS + "\n");
								bw_aln.write(qstr + "\n");
								bw_aln.write(frstr + "\n");
								bw_aln.write(refSeq + "\n");
								seqs[0] = frstr;
								seqs[1] = qstr;
								seqs[2] = refSeq;
							} else {
								String fqstr = getFlankSequence(sequence, qstart - 1, flank_size);
								String frstr = getFlankSequence(refSequence, sstart - 1, flank_size);
								if(!isPositiveStrand)
									frstr = SequenceUtil.reverseComplement(frstr); 
								int gap = qstart - 1;
								String gapSeq = sequence.substring(0, gap);
								String qstr = qrySeq.substring(0, qrySeq.length() / 2 + 1) + 
										gapSeq + 
										fqstr.substring(fqstr.length() / 2);
								bw_aln.write("@" + read_id + " AP:Z:" + (qe - fqstr.length() / 2) + "," + qe + "," + (qe + qstart) + "," + (qe + qstart + flank_size) + " RL:i:" + rlen + "\n");
								bw_aln.write("#AP:Z:" + tagAP + " FR:Z:" + tagFR + " FQ:Z:" + tagFQ + " RS:Z:" + tagRS + "\n");
								bw_aln.write("#AP:Z:" + (qe + qstart) + "," + (qe + qend) + "," + sstart + "," + send + "," + (isPositiveStrand ? "f" : "r") + 
										" FR:Z:" + frstr + "," + send + "," + (isPositiveStrand ? "f" : "r") + 
										" FQ:Z:" + fqstr + "," + qend + ",3'" + 
										" RS:Z:" + refId + "\n");
								bw_aln.write(qstr + "\n");
								bw_aln.write(refSeq + "\n");
								bw_aln.write(frstr + "\n");
								seqs[0] = refSeq;
								seqs[1] = qstr;
								seqs[2] = frstr;
							} 
							for(int i=0; i<seqs.length; i++) 
								seqs[i] = seqs[i].replaceAll("U", "T");
							
							bw_aln.write(getMSA(new String[] { seqs[0], seqs[2] }));
							bw_aln.write(getMSA(new String[] { seqs[0], seqs[1] }));
							bw_aln.write(getMSA(new String[] { seqs[2], seqs[1] }));
							bw_aln.write(getMSA(seqs));
							bw_aln.write("\n");
						} 
						if(record != null)
							read_id = record.getReadName(); 
					} 
					iterator.close();
					samReader.close();
				} 
			} 
			if(bw_complete != null) bw_complete.close(); 
			if(bw_leftover != null) bw_leftover.close(); 
			if(bw_clips != null) bw_clips.close(); 
			if(bw_qual != null) bw_qual.close(); 
			if(bw_aln != null) bw_aln.close(); 
			if(dict != null) dict.close(); 
			if(samWriter != null) samWriter.close(); 
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}

	private static String getMSA(String[] seqs) {
		// TODO Auto-generated method stub
		List<DNASequence> dnaSeqs = new ArrayList<>();
		try {
			for(String seq :seqs) {
				dnaSeqs.add(new DNASequence(seq));
			} 
		} catch (CompoundNotFoundException e) {
			e.printStackTrace();
		} 
		List<AlignedSequence<DNASequence, NucleotideCompound>> aseqs = null;
		if(seqs.length == 2) {
			SequencePair<DNASequence, NucleotideCompound> seqPair = 
					Alignments.getPairwiseAligner(dnaSeqs.get(0), dnaSeqs.get(1), 
							Alignments.PairwiseSequenceAlignerType.GLOBAL, 
							new SimpleGapPenalty(8, 4), 
							SubstitutionMatrixHelper.getNuc4_4()).getPair();
			aseqs = seqPair.getAlignedSequences();
		} else if(seqs.length > 2) {
			Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(dnaSeqs, new Object[0]);
			aseqs = profile.getAlignedSequences();
			ConcurrencyTools.shutdown();
		} 
		List<String> alignedSeqs = new ArrayList<>();
		for(AlignedSequence<DNASequence, NucleotideCompound> seq : aseqs)
			alignedSeqs.add(seq.toString()); 
		return getAlignmentString(alignedSeqs);
	}

	private static String getAlignmentString(List<String> alignedSeqs) {
		// TODO Auto-generated method stub
		StringBuilder builder = new StringBuilder();
		builder.append(alignedSeqs.get(0));
		builder.append("\n");
		for(int i = 1; i < alignedSeqs.size(); i++) {
			String a = alignedSeqs.get(i - 1);
			String b = alignedSeqs.get(i);
			for(int j = 0; j < a.length(); j++) {
				char c1 = a.charAt(j);
				char c2 = b.charAt(j);
				builder.append((c1 == c2 && c1 != '-') ? "|" : " ");
			} 
			builder.append("\n");
			builder.append(b);
			builder.append("\n");
		} 
		return builder.toString();
	}

	private static void reverseArray(byte[] arr) {
		// TODO Auto-generated method stub
		for(int left = 0, right = arr.length - 1; left < right; left++, right--) {
			byte aux = arr[left];
			arr[left] = arr[right];
			arr[right] = aux;
		} 
	}

	private static class Pair {
		private final byte k;
		private final long v;

		public Pair(byte k, long v) {
			this.k = k;
			this.v = v;
		}
	}

	private static double[] calcQSMeanAndPercentiles(String bam_file) {
		// TODO Auto-generated method stub
		SamReader samReader = SamReaderFactory.make()
				.validationStringency(ValidationStringency.LENIENT)
				.samRecordFactory((SAMRecordFactory)DefaultSAMRecordFactory.getInstance())
				.open(new File(bam_file));
		SAMRecordIterator iterator = samReader.iterator();
		Map<Byte, Long> qsCount = new HashMap<>();
		for(byte b = 0; b < 95; ) {
			qsCount.put(b, 0L);
			b = (byte)(b + 1);
		} 
		while(iterator.hasNext()) {
			SAMRecord record = (SAMRecord)iterator.next();
			if(record.getReadUnmappedFlag() || 
					record.isSecondaryOrSupplementary())
				continue; 
			int qstart = record.getReadPositionAtReferencePosition(record.getAlignmentStart());
			int qend = record.getReadPositionAtReferencePosition(record.getAlignmentEnd());
			byte[] quals = record.getBaseQualities();
			for(int i = qstart - 1; i < qend; ) {
				qsCount.put(quals[i], qsCount.get(quals[i]) + 1L);
				i++;
			} 
		} 
		iterator.close();
		try {
			samReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		} 
		long bs = 0L, qs = 0L, maxc = 0L;
		List<Pair> pairs = new ArrayList<>();
		for(byte b : qsCount.keySet()) {
			long cnt = qsCount.get(b);
			bs += cnt;
			qs += b * cnt;
			if(cnt > maxc)
				maxc = cnt; 
			pairs.add(new Pair(b, cnt));
		} 
		Collections.sort(pairs, new Comparator<Pair>() {
			public int compare(ViralChimericReadsAnalysisCmd.Pair o1, ViralChimericReadsAnalysisCmd.Pair o2) {
				return Integer.compare(o1.k, o2.k);
			}
		});
		double[] stats = new double[102];
		stats[0] = qs / bs;
		long cumBs = 0L;
		int percent = 0;
		double k1 = Double.NaN;
		for(Pair pair : pairs) {
			if(pair.v == 0L)
				continue; 
			cumBs += pair.v;
			double k = pair.k;
			double percentile;
			while((percentile = percent * 0.01D * (bs - 1L)) <= (cumBs - 1L)) {
				if(percentile < (cumBs - pair.v)) {
					long n1 = cumBs - pair.v;
					long n = bs - 1L;
					double p1 = (double)(n1 - 1L) / n;
					double p2 = (double)n1 / n;
					stats[percent + 1] = k1 + (k - k1) * (percent * 0.01D - p1) / (p2 - p1);
				} else {
					stats[percent + 1] = k;
				} 
				percent++;
			} 
			k1 = k;
		} 
		System.err.println("#total aligned bases: " + bs);
		System.err.println("#QS average: " + stats[0]);
		System.err.println("#QS quartiles: ");
		System.err.println("     0%: " + stats[1]);
		System.err.println("    25%: " + stats[26]);
		System.err.println("    50%: " + stats[51]);
		System.err.println("    75%: " + stats[76]);
		System.err.println("    99%: " + stats[100]);
		System.err.println("   100%: " + stats[101]);
		System.err.println("#QS histogram: ");
		for(Pair pair : pairs) {
			int n = (int)((double) pair.v / maxc * 100.0D);
			System.err.println(String.format("%2d", pair.k) + " " + 
					String.format("%12d", pair.v) + " " + 
					new String(new char[n]).replace("\000", "|"));
		} 
		return stats;
	}

	private static boolean isPolyATail(String seqstr) {
		// TODO Auto-generated method stub
		int[] b = new int[3];
		for(char c : seqstr.toCharArray()) {
			switch (Character.toUpperCase(c)) {
			case 'A':
				b[0] = b[0] + 1;
				break;
			case 'T':
				b[1] = b[1] + 1;
				break;
			case 'U':
				b[2] = b[2] + 1;
				break;
			}
		} 
		double max = seqstr.length() * 0.8D;
		return !(b[0] < max && b[1] < max && b[2] < max);
	}

	private static String getFlankSequence(String sequence, int p, int flank) {
		// TODO Auto-generated method stub
		return getFlankSequence(sequence, p, flank, 'N');
	}

	private static String getFlankSequence(String sequence, int p, int flank, char pad) {
		// TODO Auto-generated method stub
		StringBuilder flank_str = new StringBuilder();
		int readLen = sequence.length();
		flank_str.setLength(0);
		if(p < flank)
			for(int k = 0; k < flank - p; k++)
				flank_str.append(pad);  
		flank_str.append(sequence.substring(Math.max(0, p - flank), Math.min(readLen, p + flank + 1)));
		if(p + flank + 1 > readLen)
			for(int k = 0; k < p + flank + 1 - readLen; k++)
				flank_str.append(pad);
		// for better visualization only
		// will remove for a release
		// flank_str.replace(flank, flank+1, "["+flank_str.charAt(flank)+"]");
		return flank_str.toString();
	}

	public static boolean containsChimericRegion(byte[] quals, int windowSize, double qualThresh, int countThresh, String statsMethod) {
		// TODO Auto-generated method stub
		int count = 0;
		double q = -1.0D;
		switch (statsMethod.toLowerCase()) {
		case "median": 
			q = 0.5D;
			break;
		case "q25":
			q = 0.25D;
			break; 
		case "q50":
			q = 0.5D;
			break;
		case "q75":
			q = 0.75D;
			break; 
		case "mean":
			ArrayDeque<Byte> deque = new ArrayDeque<>();
			double sum = 0.0D;
			for(int i = 0; i <= quals.length; i++) {
				if(i >= windowSize) {
					double stat = sum / windowSize;
					count++;
					if(stat < qualThresh && count >= countThresh)
						return true; 
					sum -= deque.pop();
				} 
				if(i < quals.length) {
					deque.offer(quals[i]);
					sum += quals[i];
				} 
			} 
			break;
		default:
			throw new RuntimeException("!!!");
		} 
		if(q >= 0.0D) {
			SlidingWindowPercentile medianWindow = new SlidingWindowPercentile(quals, q);
			for(int i = 0; i <= quals.length; i++) {
				if(i >= windowSize) {
					double stat = medianWindow.getPercentile();
					count++;
					if(stat < qualThresh && count >= countThresh)
						return true; 
					medianWindow.remove(i - windowSize);
				} 
				if(i < quals.length)
					medianWindow.add(i); 
			} 
		} 
		return false;
	}

	public static double[] slidingWindowMeans(byte[] nums, int k) {
		// TODO Auto-generated method stub
		double[] means = new double[nums.length - k + 1];
		ArrayDeque<Byte> deque = new ArrayDeque<>();
		double sum = 0.0D;
		for(int i = 0; i <= nums.length; i++) {
			if(i >= k) {
				means[i - k] = sum / k;
				sum -= deque.pop();
			} 
			if(i < nums.length) {
				deque.offer(nums[i]);
				sum += nums[i];
			} 
		} 
		return means;
	}

	public static double[] slidingWindowPercentiles(byte[] nums, int k, double q) {
		// TODO Auto-generated method stub
		double[] percentiles = new double[nums.length - k + 1];
		SlidingWindowPercentile percentileWindow = new SlidingWindowPercentile(nums, q);
		for(int i = 0; i <= nums.length; i++) {
			if(i >= k) {
				percentiles[i - k] = percentileWindow.getPercentile();
				percentileWindow.remove(i - k);
			} 
			if(i < nums.length)
				percentileWindow.add(i); 
		} 
		return percentiles;
	}

	public static double percentile(byte[] nums, int startPos, int len, double q) {
		// TODO Auto-generated method stub
		if(startPos + len > nums.length)
			throw new RuntimeException("!!!"); 
		SlidingWindowPercentile percentileWindow = new SlidingWindowPercentile(nums, q);
		for(int i = 0; i < len; i++)
			percentileWindow.add(startPos + i); 
		return percentileWindow.getPercentile();
	}

	public static double percentile(byte[] nums, double q) {
		// TODO Auto-generated method stub
		return percentile(nums, 0, nums.length, q);
	}

	public static double median(byte[] nums, int startPos, int len) {
		// TODO Auto-generated method stub
		return percentile(nums, startPos, len, 0.5D);
	}

	public static double median(byte[] nums, double q) {
		// TODO Auto-generated method stub
		return percentile(nums, 0, nums.length, 0.5D);
	}

	public static double[] slidingWindowMedians(byte[] nums, int k) {
		return slidingWindowPercentiles(nums, k, 0.5D);
	}

	private static class SlidingWindowPercentile {
		private TreeSet<Integer> lower;
		private TreeSet<Integer> upper;
		private byte[] elements;
		private double percentile;

		public SlidingWindowPercentile(byte[] nums, double p) {
			if(p < 0.0D || p > 1.0D)
				throw new RuntimeException("!!!"); 
			elements = nums;
			Comparator<Integer> c = new Comparator<Integer>() {
				public int compare(Integer a, Integer b) {
					return elements[a] == elements[b] ? a - b : elements[a] < elements[b] ? -1 : 1;
				}
			};
			lower = new TreeSet<>(c);
			upper = new TreeSet<>(c);
			percentile = p;
		}

		public void add(int i) {
			if(lower.size() * (1.0D - percentile) <= upper.size() * percentile) {
				upper.add(i);
				int m = upper.first();
				upper.remove(m);
				lower.add(m);
			} else {
				lower.add(i);
				int m = lower.last();
				lower.remove(m);
				upper.add(m);
			} 
		}

		public boolean remove(int i) {
			return lower.remove(i) || upper.remove(i);
		}

		public double getPercentile() {
			if(lower.isEmpty() && upper.isEmpty())
				return Double.NaN; 
			if(lower.isEmpty())
				return elements[upper.first()]; 
			if(upper.isEmpty())
				return elements[lower.last()]; 
			int n1 = lower.size();
			int n2 = upper.size();
			int n = n1 + n2 - 1;
			double p1 = (double)(n1 - 1) / n;
			double p2 = (double)n1 / n;
			return elements[lower.last()]+(elements[upper.first()]-elements[lower.last()])*(percentile-p1)/(p2-p1);
		}
	}
}
