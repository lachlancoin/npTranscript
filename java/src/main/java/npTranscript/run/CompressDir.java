package npTranscript.run;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import npTranscript.cluster.Outputs;



public class CompressDir {
	/** Class used to write compressed file 
	* allows writing direct to zip, or to file, which is appeneded at end
	 *
	public static void main(String[] args){
		try{
			//if(true) System.exit(0);
			 File dir1 = new File(System.getProperties().getProperty("user.dir"));
	        	CompressDir.compress(dir1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	} */
	
	public static boolean doTrim = false;
	
	 FileOutputStream dest;
	    CheckedOutputStream checksum;
	    ZipOutputStream outS;
	    OutputStreamWriter osw;
	    SequenceOutputStream1 osw_s;
	    
	   public  File inDir;
	    int len;
	    
	    public CompressDir(File f, boolean includeFileLength) throws IOException{
	    	this.includeFileLength = includeFileLength;

	    	this.inDir = f;
	    	inDir.mkdir();
	    	len = inDir.getAbsolutePath().length()+1;
	    	 dest = new FileOutputStream(new File(inDir.getParentFile(), inDir.getName()+".zip"));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
	    }
	    
	    public CompressDir(File parent, String name,boolean includeFileLength) throws IOException{
	    	 dest = new FileOutputStream(new File(parent, name));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
	         this.includeFileLength = includeFileLength;
	    }
	    
	    
	    public void run(int  min_lines) throws IOException{
	    	if(this.currentStreams.size()>0){
	    		for(Iterator<Map.Entry<String, SequenceOutputStream1>> it = currentStreams.entrySet().iterator(); it.hasNext();){
	    			Map.Entry<String, SequenceOutputStream1> so = it.next();
	    			so.getValue().printAll();
	    		//	so.getValue().close();
	    		}
	    		currentStreams.clear();
	    	}
	    	Outputs.waitOnThreads(100);
	    	
	    	
	    	System.err.println("all done");
	    		  	//}
	    	File[] f = inDir.listFiles();
	    	int min_seqs = (int) Math.floor((double)min_lines/2.0);
	    
		    	//int min_seqs = (int) Math.floor((double)min_lines/2.0);
		    	for(int i=0; i<f.length; i++){
		    		if(!f[i].getName().startsWith(".")){
		    			if(!f[i].getName().startsWith("annot_") && doTrim){
		    				try{
		    				f[i] = SequenceOutputStream1.trim(f[i],min_seqs, false);
		    				}catch(Exception exc){
		    					System.err.println("problem with "+f[i]);
		    					exc.printStackTrace();
		    				}
		    			}
		    			
		    		}
		    	}
		    	
		    	
	    	writeAllFiles(f);
	    	
	    }
	   public void writeAll() {
	    	this.writeAllFiles(inDir.listFiles());
	    }
	   
	    private void writeAllFiles(File[] f) {
	    	try{
		    	//
		    	for(int i=0; i<f.length; i++){
		    		if(!f[i].getName().startsWith(".")){
		    			if(f[i]!=null && f[i].exists()){
		    				this.writeHeader(f[i]);
		    				this.delete(f[i]);
		    			}
		    		}
		    	}
		    	this.delete(inDir);
		    	this.close();
		    	}catch(Exception exc){
		    		System.err.println("problem with "+inDir);
		    		exc.printStackTrace();
		    	}
			
		}

		private int getLength(File file)  throws IOException{
	    	BufferedReader reader = new BufferedReader(new FileReader(file));
	    	int lines = 0;
	    	while (reader.readLine() != null) lines++;
	    	reader.close();
	    	return lines;
		}

		public void close() throws IOException{
	    	this.outS.close();
	    }
	    public void delete(File f) throws Exception{
	    	if(f.isDirectory()){
	    		File[] f1 = f.listFiles();
	    		for(int i=0; i<f1.length; i++){
	    			delete(f1[i]);
	    		}
	    	}
	    	f.delete();
	    }
	    final boolean includeFileLength;
	    public void writeHeader(File f) throws IOException{
	    	
	    		this.writeHeader(f,f.getAbsolutePath().substring(len));
	    	
	    }
	    public void writeHeader(File f, String newname) throws IOException{
	    	if(f.isDirectory()){
	    		File[] f1 = f.listFiles();
	    		for(int i=0; i<f1.length; i++){
	    			writeHeader(f1[i]);
	    		}
	    	}
	    	else{
	    		OutputStreamWriter osw1 = this.getWriter(newname, true);
	    		BufferedReader br = new BufferedReader(new FileReader(f));
	            String str = "";
	            while((str = br.readLine())!=null){
	        	  osw.write(str);osw.write("\n");
	           }
	          this.closeWriter(osw1);
	           br.close();
	    	}
	    }
	    
	  /*  public void closeWriter(SequenceOutputStream1 osw_s) throws IOException{
	    	if(osw_s==this.osw_s){
	    	   osw_s.flush();
	           outS.closeEntry();
	    	}else{
	    		boolean closed = osw_s.close();
	    		if(closed) this.currentStreams.remove(osw_s.entryname());
	    	}
	    }*/
	    
	    public void closeWriter(OutputStreamWriter osw) throws IOException{
	    	if(osw==this.osw){
	    	   osw.flush();
	           outS.closeEntry();
	    	}else{
	    		osw.close();
	    	}
	    }
	    
	    Map<String, SequenceOutputStream1> currentStreams = new HashMap<String, SequenceOutputStream1>();
	    
	    public SequenceOutputStream1 getSeqStream(String entry) throws IOException{
	    	//boolean writeDirectToZip = false;
	    	/*if(writeDirectToZip){
	    	//	if(append) throw new
		    	ZipEntry headings = new ZipEntry(entry);
			    outS.putNextEntry(headings);
			    return new SequenceOutputStream1(outS);
		    }else{*/
	    	SequenceOutputStream1 so1 = currentStreams.get(entry);
		    	if(so1==null){
		    		File f = new File(inDir, entry);
		    		 so1 =  new SequenceOutputStream1(f);
		    		currentStreams.put(so1.entryname(), so1);
		    		
		    	}
		    	return so1;
	    }
	    
	    public OutputStreamWriter getWriter(String entry, boolean writeDirectToZip)throws IOException{
	    		return getWriter(entry,writeDirectToZip, false );
	    }
	    public OutputStreamWriter getWriter(String entry, boolean writeDirectToZip, boolean append)throws IOException{
	    	if(writeDirectToZip){
	    	ZipEntry headings = new ZipEntry(entry);
		    outS.putNextEntry(headings);
		    return osw;
	    	}else{
	    		File f = new File(inDir, entry);
	    		if(append && ! f.exists()) append=false;
	    		return new OutputStreamWriter((new FileOutputStream(f,append)));
	    	}
		        
	    }
	    

	   /* public static void compress(File dir) {
			try{
						(new CompressDir(dir)).run();
			}catch(Exception exc){
				exc.printStackTrace();
			}
		} */
}
