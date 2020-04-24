package npTranscript.cluster;

import java.util.ArrayList;
/**
 * @author Lachlan Coin
 *
 */

public class CigarHash {
	
	
	
	public static boolean cluster_by_annotation = true;
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#toString()
	 */
	@Override
	public String toString(){
		return this.secondKey;
	}
	
	
	
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o){
		 return this.secondKey.equals(((CigarHash)o).secondKey);
		
	}


	
	/* (non-Javadoc)
	 * @see npTranscript.cluster.CigHash#clear()
	 */
	
	 public void clear() {
		secondKey = "";
		 }

	@Override
	public int hashCode(){
		 return secondKey.hashCode();
		
	}
	
	

	String secondKey="";
	public void setSecondKey(String string) {
		// TODO Auto-generated method stub
		this.secondKey = string;
	}
	}