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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;
import japsa.seq.ZipGFF;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;
import npTranscript.cluster.Annotation;
import npTranscript.cluster.CigarCluster;
import npTranscript.cluster.GFFAnnotation;
import npTranscript.cluster.Outputs.HDFObj;
import npTranscript.cluster.Outputs1;

/**
 * @author Lachlan Coin
 *
 */

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Analysis of coronavirus sequence data")
public class PrintResults extends CommandLine {

	public PrintResults() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());
		addString("h5file", null, "Input h5", true);
		addString("annotation", null, "ORF annotation file or GFF file", false);
		addString("optsFile", null, "Name of file with extra options", false);
		addString("resdir", "./", "results directory");//"results"+System.currentTimeMillis()
		addString("reference", null, "ref genome");
		
		addBoolean("firstIsTranscriptome", false, "whether first bam is transcriptome");
		addString("gffThresh","1", "reports if greater than in dataset");
		addInt("maxTranscriptsPerGeneInGFF",1,"Maximum number of transcripts per gene (highest abundance first");
		addString("annotType", null, "Type of annotation (only included if annotation is GFF file", false);
		addString("GFF_features", "gene_name:description:gene_id:gene_biotype:gene_id", "GFF feature names");
		addBoolean("writeGFF", true, "whether to output gff ");

		
	}
	public static CigarCluster readH5(IHDF5Reader reader, String chrom, String key){
			String id1 = "/transcripts/"+chrom+"/"+key;
			
			List<String > str1 = reader.getGroupMembers(id1);
			HDFObj[] objs = new HDFObj[str1.size()];
			for(int j=0; j<str1.size(); j++){
				String id2 = id1+"/"+str1.get(j);
				objs[j] =  new  HDFObj(reader.readString(id2));
				
			}
			return new CigarCluster(chrom,key, str1.toArray(new String[0]),objs);
			//System.err.println("h");
	}
	static class ClusterIterator implements Iterator<CigarCluster> {
		final IHDF5Reader reader;
		String chrom;
		final Iterator<String> chroms;
		Iterator<String> keys;
		final String[] sources;
		ClusterIterator(String h5file){
			 reader= HDF5Factory.openForReading(h5file);
			sources =  reader.readStringArray("header");
			chroms = new HashSet<String>(reader.getGroupMembers("/transcripts")).iterator();
			chrom = chroms.next();
			keys= new HashSet<String>(reader.getGroupMembers("/transcripts/"+chrom)).iterator();
		}
		public void close(){
			reader.close();
		}
		
			@Override
			public boolean hasNext() {
				return keys.hasNext() || chroms.hasNext();
			}
			@Override
			public CigarCluster next() {
				if(!keys.hasNext()){
					chrom = chroms.next();
					keys= new HashSet<String>(reader.getGroupMembers("/transcripts/"+chrom)).iterator();
				}
				return readH5(reader, chrom, keys.next());
			}
	}
	public static void main(String[] args1) throws IOException, InterruptedException {
		CommandLine cmdLine = new PrintResults();
		String[] args2 = ViralTranscriptAnalysisCmd2.readOpts(cmdLine, args1,"optsFile");
		
		Outputs1.firstIsTranscriptome = cmdLine.getBooleanVal("firstIsTranscriptome");
		String[] gffThresh_1 = cmdLine.getStringVal("gffThresh").split(":");
		Outputs1.gffThresh = new int[gffThresh_1.length];
		
		//Outputs.gffThreshTranscriptSum=0;
		
			for(int k=0; k<gffThresh_1.length; k++){
				Outputs1.gffThresh[k] = Integer.parseInt(gffThresh_1[k]);
			}
		
		Outputs1.maxTranscriptsPerGeneInGFF = cmdLine.getIntVal("maxTranscriptsPerGeneInGFF");
		Outputs1.writeGFF = cmdLine.getBooleanVal("writeGFF");

		
		String resdir = cmdLine.getStringVal("resdir");
		File resDir = new File(resdir);
		if(!resDir.exists())resDir.mkdir();
		
		
		
		String h5file = cmdLine.getStringVal("h5File");
		ClusterIterator it = new ClusterIterator(h5file);
		
		
		String annot_file = cmdLine.getStringVal("annotation");
		int len = it.sources.length;

		String reference = cmdLine.getStringVal("reference");
		
		
		if(false){
			File annotSummary = new File(resdir, "annotation.csv.gz");
			if(annotSummary.exists()) annotSummary.delete();
			PrintWriter annotation_pw = new PrintWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(annotSummary, false))));
			File gffFile=annot_file==null ? null : new File(annot_file);
			Annotation annot  = null;
			ZipFile anno = null;
			if(gffFile!=null && (gffFile.getName().indexOf(".gff")>=0 || gffFile.getName().indexOf(".gtf")>=0)){
					anno  = getAnnotfile(gffFile);
					annot = new GFFAnnotation(anno, annotation_pw, gffFile.getName().indexOf(".gff")<0);
			}else{
				annot = 
					new Annotation(new File(annot_file),  annotation_pw, len);
			}
			if(anno!=null) anno.close();
			if(annotation_pw!=null) annotation_pw.close();
		}
		//String refFile = cmdLine.getStringVal("reference");
		ArrayList<Sequence> genomes = reference==null ? null : SequenceReader.readAll(reference, Alphabet.DNA());
		Map<String, Integer> m = new HashMap<String, Integer>();
		if(genomes!=null){
			for(int i=0; i<genomes.size(); i++){
				m.put(genomes.get(i).getName(), i);
			}
		}
	//	altT.
		
		//o.writeTotalString(this);
		
		/*	while(  it.hasNext()) {
					CigarCluster nxt = it.next();
					process1(nxt, o);//, chrom, chrom_index, geneNames);
					int  totalDepth = nxt.readCountSum();
				//	if(Outputs.writeIsoforms) o.writeIsoforms(nxt, this, totalDepth);
					o.writeDepthH5(nxt, this,  totalDepth);
			}*/
		Outputs1 	outp = new Outputs1(resDir,  it.sources); 
		while(it.hasNext()){
			CigarCluster cc = it.next();
			Sequence genome = genomes==null ? null :genomes.get(m.get(cc.chrom));
			cc.process1(outp, genome);
			System.err.println(cc.chrom);
		}
		outp.close();
		it.close();
	}
	private static ZipFile getAnnotfile(File gffFile) throws ZipException, IOException {
		ZipFile anno = null;
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
		return anno;
	}
}
