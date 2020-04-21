package npTranscript.run;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import japsa.util.CommandLine;
import japsa.util.deploy.Deployable;

@Deployable(scriptName = "npTranscript.run", scriptDesc = "Clustering of results from ViralTranscriptAnalysisCmd2")
public class ExtractClusterCmd extends CommandLine {
	
	
	public ExtractClusterCmd() {
		super();
		Deployable annotation = getClass().getAnnotation(Deployable.class);
		setUsage(annotation.scriptName() + " [options]");
		setDesc(annotation.scriptDesc());

		addString("inDir", null, "Name of inputDir", true);
		addStdHelp();
		
	}
	public void readTranscripts(File transcriptsF) throws IOException{
		transcripts = new Table(transcriptsF,"\t");
	}
	
	
	
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new ExtractClusterCmd();
		args = cmdLine.stdParseLine(args);

		String inDir = cmdLine.getStringVal("inDir");
		File in = new File(inDir);
		File transcriptF = new File(in, "0transcripts.txt.gz");
		
		File readsF = new File(in, "0readToCluster.txt.gz");
		String chrom = "0";
		
		try{
			ExtractClusterCmd ec = new ExtractClusterCmd();
			ec.readTranscripts(transcriptF);	
			String[] types = ec.transcripts.getColElements("type_nme");
			String sort = "countTotal";
		//	String[] filters = "type_nme;5_3".split(":");
			String[] filters1 ="countTotal;10;true".split(":"); 
			String[] grpName = "upstream:downstream".split(":");
			for(int i=0; i<types.length; i++){
				String[] filters = new String[] {"type_nme",types[i]};
				CompressDir outdir = new CompressDir(in, "clusters_"+types[i]+".zip");
				ec.run(readsF, outdir,	sort, filters, 
					filters1,grpName, false);
			}
			
			
			
			
			//ec.run(transcriptF, readsF, outdir,	"countTotal", "type_nme;5_3".split(":"), 
				//	"countTotal;10;true".split(":"),"upstream:downstream".split(":"));
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
			if(col<0) throw new RuntimeException("!!");
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
	Table transcripts;
	CompressDir outdir;
	
	
	//
	public void run(  File readsF,  CompressDir outdir,
			String sort, String[] filters, String[] limits,String[] groupIds, boolean print
			) throws IOException{
		//if(!outdir.exists()) outdir.mkdir();
		reads = filters==null ? new Table(readsF,"\t") : new Table(readsF,"\t", filters);
		Table transcripts1 = new Table(this.transcripts);
	    this.outdir = outdir;
		if(sort!=null) transcripts1.sort(sort, true);
		
		if(limits!=null){
			for(int i=0; i<limits.length; i++){
				String[] filt = limits[i].split(";");
				transcripts1.limit(filt[0], Integer.parseInt(filt[1]), Boolean.parseBoolean(filt[2]));
			}
		}
		String[] ids = transcripts1.getCol("ID");
		Map<String, String> m= new HashMap<String, String>();
		for(int i=0; i<ids.length; i++){
			List<String> all_ids = Arrays.asList(ids[i].split(";"));
			Table t2 = this.reads.extract("clusterId",groupIds, all_ids);
			if(t2.data.size()>0){
			t2.print(outdir);
			for(Iterator<String> it = all_ids.iterator(); it.hasNext();){
				m.put(it.next(), t2.name);
			}
			}
//			t2.print(outdir);
		}
		transcripts1.replaceAll("ID", m);
		if(print) transcripts1.print(outdir);
		
	}
	public static void cluster(String resdir,  File transcriptF, File readsF, 
			boolean filterByType) {
		File in = new File(resdir);
		//File transcriptF;// = profiles.get(i).transcripts_file;
		//File readsF = ;
		try{
			ExtractClusterCmd ec = new ExtractClusterCmd();
			ec.readTranscripts(transcriptF);	
			String[] types = ec.transcripts.getColElements("type_nme");
			String sort = "countTotal";
			String[] filters1 ="countTotal;1;true".split(":"); 
			String[] grpName = "chrom:upstream:downstream".split(":");
			if(filterByType){
			for(int ij=0; ij<types.length; ij++){
				String[] filters = new String[] {"type_nme",types[ij]};
				File outdir = new File(in, "clusters_"+types[ij]);
				CompressDir comp = new CompressDir(in, "clusters_"+types[ij]+".zip");
				ec.run(readsF, comp,	sort, filters, 
					filters1,grpName, false);
				comp.close();
				//CompressDir.compress(outdir);
			}
			}else{
			//	File outdir = new File(in, "clusters");
				CompressDir comp = new CompressDir(in, "clusters.zip");
				ec.run(readsF, comp,	sort, null, filters1,grpName, false);
				comp.close();
			}
		}catch(Exception exc){
			exc.printStackTrace();
		}
		
	}
	

		
		
	}


