package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Filter fastq")

public class FilterFastq extends CommandLine {
	 static final FastqWriterFactory factory = new FastqWriterFactory();

	 public FilterFastq(){
			super();
			Deployable annotation = getClass().getAnnotation(Deployable.class);		
			setUsage(annotation.scriptName() + " [options]");
			setDesc(annotation.scriptDesc());

			addString("fastq", null,  "Name of fastq file", true);
			addString("readList", null,  "reads to include", false);
	 }
	 static Collection[] getReads(String[] readList) throws IOException {
			Collection[] reads = null;
			if(readList!=null && readList.length>0){
			reads = new Collection[readList.length];
			 for(int i=0; i<reads.length; i++){
				reads[i] = new HashSet<String>();
				BufferedReader br = new BufferedReader(new FileReader(new File(readList[i])));
				String st;
				while((st = br.readLine())!=null){
					String st_ = st.split("\\s+")[0];
				//	System.err.println(st_);
					reads[i].add(st_);
				}
				br.close();
			 }
			}
			return reads;
			}
	 public static void main(String [] args) throws IOException, InterruptedException{		 		
			CommandLine cmdLine = new FilterFastq();		
			args = cmdLine.stdParseLine(args);		

			String bamFile = cmdLine.getStringVal("fastq");
			String readL= cmdLine.getStringVal("readList");
			
			File resDir = new File("./results");
			resDir.mkdir();
			String[] bamFiles_ = bamFile.split(":");
			if(bamFile.equals("all") || bamFile.equals(".")){
				bamFiles_ = (new File("./")).list(new FilenameFilter(){

					@Override
					public boolean accept(File dir, String name) {
						return name.endsWith(".fastq");
					}
					
				});
			}
			Collection[]reads = getReads(readL.split(":"));
			FastqWriter fastq = factory.newWriter( new File(resDir,"combined.fastq"));
			for(int i=0; i<bamFiles_.length; i++){
				fastqFilter( new File(bamFiles_[i]), reads, fastq);		
			}
			fastq.close();

			//paramEst(bamFile, reference, qual);
		}

	 
	
	 
	public static void fastqFilter( File input, Collection[] reads, FastqWriter fastq){
	//	if(resDir.equals(input.getParentFile())) throw new RuntimeException("output and input dir same!!");
		FastqReader reader = (new FastqReader(input));
		
		//File out = new File(resDir,input.getName());
		//System.err.println(out.getAbsolutePath());
		//if(out.exists()) throw new RuntimeException("!!");
	
		outer: while (reader.hasNext()){
			FastqRecord seq = reader.next();
			for(int i=0; i<reads.length; i++){
				String rn = seq.getReadName().split("\\s+")[0];
			if(reads[i].contains(rn)){
				fastq.write(seq);
				continue outer;
			}
			}
		}
		reader.close();
	}
	
	
}
