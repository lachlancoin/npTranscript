package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Merging results from ViralTranscriptAnalysisCmd2")
public class MergeCmd extends CommandLine {
	
	
	public MergeCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("inDir","./", "Name of inputDir", false);
		addString("tocombine", "", "Directories to merge", false);
		//"\t5_3\t"
		addString("filter", null, "type to filter, e.g. \t5_3\t");
		addString("ignore",null, "pattern to ignore", false);
		addString("resdir","./merged", "resultsdirectory", false);
		addStdHelp();
		
	}
	
	//ID	start	end	type_nme	breaks	hash	startBreak	endBreak	leftGene	rightGene	totLen	countTotal	count0count1	depth0	depth1	errors0	errors1	error_ratio0	error_ratio1
//  	9	8	29858	5_3	0,1,2624,2989	-1287377339	69	26237	LEADER	E	0	3	3	0	11012	0	707	00.0642	NaN

	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		MergeCmd cmdLine = new MergeCmd();
		args = cmdLine.stdParseLine(args);
		String tf_name = "0transcripts.txt.gz";
		String pattern = cmdLine.getStringVal("ignore");
		String filter = cmdLine.getStringVal("filter");
		String resdir = cmdLine.getStringVal("resdir");
		//args.cmdLine.setOutDir(resdir);
		File outdir = new File(resdir);
		if(!outdir.exists()) outdir.mkdir();
		
		final List<String> tocombine = Arrays.asList(cmdLine.getStringVal("tocombine").split(":"));
		
		File inDir_ = new File(cmdLine.getStringVal("inDir"));
		List<File> inDir = Arrays.asList(inDir_.listFiles(new FileFilter(){

			@Override
			public boolean accept(File pathname) {
				if(tocombine.size()>1){
					return pathname.isDirectory() && tocombine.contains(pathname.getName());
				}
				else return pathname.isDirectory() 
						&& (new File(pathname, tf_name).exists() 
						&& !pathname.getName().contains("merged")
						&& ( pattern==null || !pathname.getName().contains(pattern)));
			}
			
		}));
		Collections.sort(inDir);
		/*File[] transcriptF = new File[inDir.length];
		for(int i=0; i<inDir.length; i++){
			transcriptF[i] = new File(inDir[i], tf_name);
		}*/
	
		
		try{
			MergeCmd ec = new MergeCmd();
			ec.outdir = outdir;
			
			ec.run(inDir.toArray(new File[0]), tf_name, filter);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
  static int cntr = 1; 
  static String split = "\t";
	static class StringComp implements Comparator<String[]>{
		final int col;
		final int mult;
		
		StringComp(int col, boolean decr){
			this.col = col;
			mult = decr ? -1 : 1;
		}
		@Override
		public int compare(String[] o1, String[] o2) {
			Integer i1 = Integer.parseInt(o1[col]);
			Integer i2 = Integer.parseInt(o2[col]);
			return mult*i1.compareTo(i2);
		}

		
		
	}
	Table reads;
	Table[] transcripts;
	 File outdir ;
	
	Map<String, String[][]> map = new HashMap<String, String[][]>();
	
	
	static class MergeInds{
		List<Integer> remaining = new ArrayList<Integer>();
		Integer[] merge_inds;
		private String[] target;
		int[] merge_offset;
		String[][] head_input;
		final String default_val;
		public String[] merge( String[][]input , boolean add_j ){
			
			boolean first= true;
			for(int j=0; j<input.length; j++){
				String[] str = input[j];
				if(str!=null && first){
					first = false;
					for(int i=0; i<remaining.size(); i++){
						target[i] = str[remaining.get(i)];
					}
				}
				for(int k=0; k<merge_inds.length; k++){
						String nme = str==null ? default_val : str[merge_inds[k]];
						if(add_j) nme = nme.replaceAll("0" , ""+j);
						target[merge_offset[k]+j] = nme;
				}
			}
			return target.clone();
		}
		public String[] newHeader(){
			return this.merge(head_input, true);
		}
		
		public MergeInds(List<String> header, String[] merge, String[] toignore, int num_sources, String default_val){
			merge_inds = new Integer[merge.length];
			
			merge_offset = new int[merge.length];
			this.default_val = default_val;
			for(int i=0; i<header.size(); i++){
				remaining.add(i);
			}
				
			for(int i=0; i<merge_inds.length; i++){
				merge_inds[i] = header.indexOf(merge[i]);
				if(merge_inds[i]<0) {
					throw new RuntimeException("could not find "+merge[i]);
				}
				remaining.remove(merge_inds[i]);
			}
			for(int i=0; i<toignore.length; i++){
				Integer torem = header.indexOf(toignore[i]);
				if(torem<0) throw new RuntimeException("!!");
				remaining.remove(torem);
			}
			target = new String[remaining.size() + num_sources*merge.length];
			{
			int i= remaining.size();
			for(int k=0; k<merge_inds.length; k++){
					merge_offset[k] = i;
					i+=num_sources;
			}
			}
			 head_input  = new String[num_sources][];
			for(int i=0; i<num_sources; i++){
				head_input[i] = header.toArray(new String[0]);
			}
		
		//	merge(merged_header, input, true);
			
		}
	}
	static void mkSymLink(File src, File target){
		if(target.exists()){
			System.err.println("exists");
		return;
		}
		try{
		
		String cmd = "ln -s "+src.getAbsolutePath()+" "+target.getAbsolutePath();
		Process process = Runtime.getRuntime().exec(cmd);
		process.waitFor();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	static void print(File in, PrintWriter out, String nme, String val, boolean first, String pattern) throws IOException{
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(in))));
		String head = br.readLine();
		if(first){
			out.println(head);
		}
		List<String> header = Arrays.asList(head.split("\t"));
		int ind = header.indexOf(nme);
		if(ind<0) {
			throw new RuntimeException("!!");
		}
		String str = "";
		
		boolean same = val.equals("0");
		while((str = br.readLine())!=null){
			if(pattern == null || str.indexOf(pattern)>=0){
				if(same) out.println(str);
				else{
					String[] st_ = str.split("\t");
					st_[ind] = val;
					out.print(st_[0]);
					for(int j=1; j<st_.length ; j++){
						out.print("\t");out.print(st_[j]);
					}
					out.println();
				}
			}
		}
		br.close();
	}
	
	public void run(File[] inDir , String tf_name, String filter) throws IOException{
		String[] ids = new String[] {"type_nme", "leftGene", "rightGene"};
		System.err.println(Arrays.asList(inDir));
		List<String> header = null;
		int num_sources = inDir.length;
		StringBuffer info = new StringBuffer(); 
		
		
		// this combines the reads file
		PrintWriter reads_pw  = new PrintWriter(
				new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outdir, "0readToCluster.txt.gz")))));
		for(int i=0; i<inDir.length; i++){	
			//File  targetFile = new File(outdir, "0readToCluster."+i+".txt.gz");
			File readsFile = new File(inDir[i], "0readToCluster.txt.gz");
			System.err.println(readsFile.getAbsolutePath());
			print(readsFile, reads_pw,"source", (i+""), i==0, filter);
			info.append(inDir[i].getName());
			 info.append("\t");
			Table transcripts = new Table(new File(inDir[i], tf_name),"\t");
			header = transcripts.header;
			transcripts.addToMap(map, ids, i, num_sources);
		}
		reads_pw.close();
		
		
		String[] merge = new String[] {"count0",  "depth0", "errors0", "error_ratio0"};
		String[] toignore = new String[] {"countTotal"};
		MergeInds mi = new MergeInds(header,  merge,  toignore, num_sources, "0");
		List<String> head_new = Arrays.asList(mi.newHeader());
		Table t1 = new Table(head_new, tf_name);
		for(Iterator<String[][] > it = map.values().iterator(); it.hasNext();){
			String[]str = mi.merge(it.next(), false);
			t1.data.add(str);
		}
		
		t1.print(outdir, info.toString());
		
	}
	

		
		
	}


