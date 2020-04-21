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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import npTranscript.run.ExtractClusterCmd.StringComp;

public class Table{
	
	
	
	
	
	List<String> header;
	List<String[]> data = new ArrayList<String[]>();
	String name;
	
	
	public Table(List<String> header, String name){
		this.header = header;
		this.name  = name; 
	}
	//this reuses the elements in data
	public Table(Table in ){
		this.name = new String(in.name);
		header =in.header;
		for(int i=0; i<in.data.size(); i++){
			this.data.add(in.data.get(i));
		}
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
	public Table(File in, String split , String[] filter) throws IOException{
		this.name = in.getName();
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(in))));
		String st = "";
		 while((st = br.readLine()).startsWith("#")){
			 
		 }
		 header = Arrays.asList(st.split(split));
		 int col = header.indexOf(filter[0]);
			
			
		 while((st = br.readLine())!=null){
			 String[] str = st.split(split);
			 if(str[col].equals(filter[1]))data.add(str);
		 }
	}
	
		//transcripts1.filter(filters[0], filters[1]);
	

	
	
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
	
	public String[] getColElements(String nme){
		int col = header.indexOf(nme);
		Set<String>s = new HashSet<String>();
		for(int i=0; i<data.size(); i++){
			s.add(data.get(i)[col]);
		}
		return s.toArray(new String[0]);
	}
	
	
	public Table  extract(String mainID, String[] grp_ids, List<String> vals ){
		//String nme, String lab1, String lab2
		String id = ExtractClusterCmd.cntr+ "";
		int main_ind = header.indexOf(mainID);
		int[] col = new int[grp_ids.length];
		for(int i=0; i<col.length; i++){
			col[i] = header.indexOf(grp_ids[i]);
			if(col[i]<0) throw new RuntimeException("did not find "+grp_ids[i]);
		}
		ExtractClusterCmd.cntr++;
		Table t1 = new Table(this.header, id);
		String label = null;
		for(int i=0; i<this.data.size();  i++){
			String[] row = data.get(i);
			String v = (row[main_ind]);
			//System.err.println(v);
			if(vals.indexOf(v)>=0 ){
				//if(i==0){
				if(label==null){
				StringBuffer sb = new StringBuffer();
				for(int ij=0; ij<col.length; ij++){
					if(col[ij]<row.length)sb.append(row[col[ij]]+"_");
					else{
						System.err.println("warning "+ col[ij]+" is longer than row length "+Arrays.asList(row));
					}
				}
				label = sb.toString();
				}
					//lab1_ = row[col1];
					//lab2_ = row[col2];
				//	}
				t1.data.add(row);
			}
			
		}
		//if(t1.data.size()==0) throw new RuntimeException("could not find anything");
		t1.name = label+t1.data.size();
		return t1;
	}
	public void print(CompressDir outdir) throws IOException{
		print(outdir, null);
	}
	public void print(CompressDir outdir, String info) throws IOException{
		String name1 = this.name;
		OutputStreamWriter transcripts_pw = outdir.getWriter(name1, true);
		if(info!=null){
			transcripts_pw.write("#");transcripts_pw.write(info);transcripts_pw.write("\n");
		}
		transcripts_pw.write(getStr(this.header.toArray(new String[0]), ExtractClusterCmd.split));
		transcripts_pw.write("\n");
		for(int i=0; i<this.data.size(); i++){
			transcripts_pw.write(getStr(data.get(i), ExtractClusterCmd.split));
			transcripts_pw.write("\n");
		}
		outdir.closeWriter(transcripts_pw);
		
	}
	public void print(File outdir, String info) throws IOException{
		if(!outdir.exists()) outdir.mkdir();
		String name1 = 	this.name.endsWith(".gz") ? this.name :  this.name+".gz";
		
		PrintWriter transcripts_pw  = new PrintWriter(
				new OutputStreamWriter(new GZIPOutputStream(
					new FileOutputStream(new File(outdir, name1)
								))));
			
		if(info!=null){
			transcripts_pw.print("#");transcripts_pw.println(info);
		}
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