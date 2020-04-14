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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import npTranscript.run.ExtractClusterCmd.StringComp;
import npTranscript.run.MergeCmd.MergeInds;

public class Table{
	
	
	
	
	
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
		String st = "";
		 while((st = br.readLine()).startsWith("#")){
			 
		 }
		 header = Arrays.asList(st.split(split));
		
		 while((st = br.readLine())!=null){
			 data.add(st.split(split));
		 }
	}
	
	
	//public void merge(List<String>)
	//countTotal
	void sort(String nme, boolean decr){
		this.name = name+".sort="+nme;
		int col = header.indexOf(nme);
		if(col<0) throw new RuntimeException("not found "+nme+ " in "+header.toString());
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
		String id = ExtractClusterCmd.cntr+ "";
		
		
		int col = header.indexOf(nme);
		int col1 = header.indexOf(lab1);
		int col2 = header.indexOf(lab2);
		//System.err.println(header);
		//System.err.println(nme+" "+col);
		ExtractClusterCmd.cntr++;
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
	public void print(File outdir) throws IOException{
		print(outdir, null);
	}
	public void print(File outdir, String info) throws IOException{
		String name1 = this.name.endsWith(".gz") ? this.name :  this.name+".gz";
	
		PrintWriter transcripts_pw  = new PrintWriter(
				new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outdir, name1)))));
		if(info!=null) transcripts_pw.print("#");transcripts_pw.println(info);
		transcripts_pw.println(getStr(this.header.toArray(new String[0]), ExtractClusterCmd.split));
		for(int i=0; i<this.data.size(); i++){
			transcripts_pw.println(getStr(data.get(i), ExtractClusterCmd.split));
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
	public void addToMap(Map<String, String[][]> map, String[] ids, int pos, int num_sources) {
		int[] inds = new int[ids.length];
		for(int i=0; i<inds.length; i++){
			inds[i] = this.header.indexOf(ids[i]);
		}
		Iterator<String[]> it = this.data.iterator();
		while(it.hasNext()){
			String[] str = it.next();
			String key = getKey(str, inds);
			String[][] l = map.get(key);
			if(l==null) map.put(key, l = new String[num_sources][]);
			l[pos] = str;
		}
		
	}
	private static String getKey(String[] str, int[] inds) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<inds.length; i++) {
			sb.append(str[inds[i]]);
			sb.append(",");
		}
		return sb.toString();
	}
	
	
}