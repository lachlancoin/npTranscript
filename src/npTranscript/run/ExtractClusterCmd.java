package npTranscript.run;

import java.io.BufferedReader;
import java.io.File;
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
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

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
	
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLine cmdLine = new ExtractClusterCmd();
		args = cmdLine.stdParseLine(args);

		String inDir = cmdLine.getStringVal("inDir");
		File in = new File(inDir);
		File transcriptF = new File(in, "0transcripts.txt.gz");
		File outdir = new File(in, "clusters");
		outdir.mkdir();
		File readsF = new File(in, "0readToCluster.txt.gz");
	
		
		try{
			ExtractClusterCmd ec = new ExtractClusterCmd();
			ec.run(transcriptF, readsF, outdir);
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
	 class Table{
		List<String> header;
		List<String[]> data = new ArrayList<String[]>();
		String name;
		public Table(List<String> header, String name){
			this.header = header;
			this.name  = name; 
		}
		public Table(File in, String split ) throws IOException{
			this.name = in.getName();
			BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(in))));
			 header = Arrays.asList(br.readLine().split(split));
			 String st = "";
			 while((st = br.readLine())!=null){
				 data.add(st.split(split));
			 }
		}
		//countTotal
		void sort(String nme, boolean decr){
			this.name = name+".sort="+nme;
			int col = header.indexOf(nme);
			StringComp sc = new StringComp(col, decr);
			Collections.sort(data, sc);
		}
		public void filter(String nme, String val){
			this.name = name+".filter="+nme;
			int col = header.indexOf(nme);
			
			for(int i=this.data.size()-1; i>=0; i--){
				if(!data.get(i)[col].equals(val)) data.remove(i);
			}
		}
		public void limit(String nme, int lim, boolean greater) {
			// TODO Auto-generated method stub
			this.name = name+".limit="+nme;
			int col = header.indexOf(nme);
			for(int i=this.data.size()-1; i>=0; i--){
				int v = Integer.parseInt(data.get(i)[col]);
				if(greater && v < lim || !greater && v > lim ){
					data.remove(i);
				}
			}
			
		}
		
		public void replaceAll(String nme, Map<String, String> m) {
			this.name = name+".replace="+nme;
			int col = header.indexOf(nme);
			for(int i=this.data.size()-1; i>=0; i--){
				String[] row = data.get(i);
				String[] key = row[col].split(";");
				String val = m.get(key[0]);
				if(val!=null) row[col] =val;
				
			}
			
		}
		
		public String[] getCol(String nme){
			int col = header.indexOf(nme);
			String[] out = new String[this.data.size()];
			for(int i=0; i<out.length; i++){
				out[i] = data.get(i)[col];
			}
			return out;
		}
		
		
		
		public Table  extract(String nme, String lab1, String lab2, List<String> vals ){
			String id = cntr+ "";
			
			
			int col = header.indexOf(nme);
			int col1 = header.indexOf(lab1);
			int col2 = header.indexOf(lab2);
			//System.err.println(header);
			//System.err.println(nme+" "+col);
			cntr++;
			Table t1 = new Table(this.header, id);
			String lab1_="" ;
			String lab2_="";
			for(int i=0; i<this.data.size();  i++){
				String[] row = data.get(i);
				String v = (row[col]);
				//System.err.println(v);
				if(vals.indexOf(v)>=0){
					//if(i==0){
						lab1_ = row[col1];
						lab2_ = row[col2];
					//	}
					t1.data.add(row);
				}
				
			}
			if(t1.data.size()==0) throw new RuntimeException("could not find anything");
			t1.name = lab1_+"_"+lab2_+"_"+t1.data.size();
			return t1;
		}
		public void print() throws IOException{
			PrintWriter transcripts_pw  = new PrintWriter(
					new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outdir.getParentFile(), this.name+".gz")))));
			transcripts_pw.println(getStr(this.header.toArray(new String[0]), split));
			for(int i=0; i<this.data.size(); i++){
				transcripts_pw.println(getStr(data.get(i), split));
			}
			transcripts_pw.close();
		}
		
		
		private String getStr(String[] header2, String split) {
			// TODO Auto-generated method stub
			StringBuffer sb = new StringBuffer();
			sb.append(header2[0]);
			for(int i=1; i<header2.length; i++){
				sb.append(split); sb.append(header2[i]);
			}
			return sb.toString();
		}
		
		
	}
	Table reads;
	Table transcripts;
	File outdir;
	
	
	
	public void run(File transcriptsF , File readsF,  File outdir) throws IOException{
		reads = new Table(readsF,"\t");
		transcripts = new Table(transcriptsF,"\t");
	    this.outdir = outdir;
		transcripts.sort("countTotal", true);
		transcripts.filter("type_nme",  "5_3");
		transcripts.limit("countTotal", 10, true);
		String[] ids = transcripts.getCol("ID");
		Map<String, String> m= new HashMap<String, String>();
		for(int i=0; i<ids.length; i++){
			List<String> all_ids = Arrays.asList(ids[i].split(";"));
			Table t2 = this.reads.extract("clusterId", "upstream", "downstream", all_ids);
			t2.print();
			for(Iterator<String> it = all_ids.iterator(); it.hasNext();){
				m.put(it.next(), t2.name);
			}
			t2.print();
		}
		transcripts.replaceAll("ID", m);
		transcripts.print();
		
	}
	

		
		
	}


