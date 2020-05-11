package npTranscript.run;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.zip.Adler32;
import java.util.zip.CheckedOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import japsa.seq.SequenceOutputStream;



public class CompressDir {
	/** Class used to write compressed file 
	* allows writing direct to zip, or to file, which is appeneded at end
	 * */
	public static void main(String[] args){
		try{
			//if(true) System.exit(0);
			 File dir1 = new File(System.getProperties().getProperty("user.dir"));
	        	CompressDir.compress(dir1);
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
	 FileOutputStream dest;
	    CheckedOutputStream checksum;
	    ZipOutputStream outS;
	    OutputStreamWriter osw;
	    SequenceOutputStream osw_s;
	    
	   public  File inDir;
	    int len;
	    
	    public CompressDir(File f) throws IOException{
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
	    
	    public CompressDir(File parent, String name) throws IOException{
	    	 dest = new FileOutputStream(new File(parent, name));
	         checksum = new   CheckedOutputStream(dest, new Adler32());
	         outS = new 
	         ZipOutputStream(new 
	           BufferedOutputStream(checksum));
	         osw = new OutputStreamWriter(outS);
	         outS.setMethod(ZipOutputStream.DEFLATED);
	    }
	    
	    
	    public void run() throws IOException{
	    	try{
	    	File[] f = inDir.listFiles();
	    	for(int i=0; i<f.length; i++){
	    	this.writeHeader(f[i]);
	    	this.delete(f[i]);
	    	}
	    	this.delete(inDir);
	    	this.close();
	    	}catch(Exception exc){
	    		System.err.println("problem with "+inDir);
	    		exc.printStackTrace();
	    	}
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
	    
	    public void closeWriter(SequenceOutputStream osw_s) throws IOException{
	    	if(osw_s==this.osw_s){
	    	   osw_s.flush();
	           outS.closeEntry();
	    	}else{
	    		osw_s.close();
	    	}
	    }
	    
	    public void closeWriter(OutputStreamWriter osw) throws IOException{
	    	if(osw==this.osw){
	    	   osw.flush();
	           outS.closeEntry();
	    	}else{
	    		osw.close();
	    	}
	    }
	    public SequenceOutputStream getSeqStream(String entry,  boolean append) throws IOException{
	    	boolean writeDirectToZip = false;
	    	if(writeDirectToZip){
	    	//	if(append) throw new
		    	ZipEntry headings = new ZipEntry(entry);
			    outS.putNextEntry(headings);
			    return new SequenceOutputStream(outS);
		    	}else{
		    		File f = new File(inDir, entry);
		    		if(append && ! f.exists()) append=false;
		    		return new SequenceOutputStream((new FileOutputStream(f, append)));
		    	}
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

	    public static void compress(File dir) {
			try{
						(new CompressDir(dir)).run();
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}
}
