package npTranscript.cluster;

/**
 * @author Lachlan Coin
 *
 */

public class CigarHash  {
	
	
	
	public static boolean cluster_by_annotation = true;
	public CigarHash(String string, int i) {
		this.secondKey = string;
		if(secondKey.length()==0){
			System.err.println("warning");
		}
		//this.end = i;
	}



	public CigarHash() {
		// TODO Auto-generated constructor stub
	}



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
	
	

	String secondKey="default";
//	Integer  end;
	public void setSecondKey(String string) {
		// TODO Auto-generated method stub
	//	this.end = end;
		
		this.secondKey = string;
		if(secondKey.length()==0){
			System.err.println("warning");
		}
	}



	/*@Override
	public int compareTo(Object o) {
		// TODO Auto-generated method stub
		if(((CigarHash)o).secondKey.equals(this.secondKey)) return 0;
		int res = end.compareTo(((CigarHash)o).end);
		if(res==0){ // we may not need this, it might be ok if two hashes are equal under compare, but to be sure
			//secondKey.com
			return Integer.compare(secondKey.hashCode(), ((CigarHash)o).hashCode());
		}
		return res;
	}*/
	}