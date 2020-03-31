package npTranscript.cluster;

import java.util.ArrayList;
/**
 * @author Lachlan Coin
 *
 */

public class CigarHash extends ArrayList<Integer>{
		public static int round = 100;
		@Override
		public String toString(){
			StringBuffer sb = new StringBuffer();
			for(int i=0; i<this.size(); i++){
				if(i>0) sb.append(",");
				sb.append(this.get(i));
			}
			return sb.toString();
		}
		
		
		@Override
		public boolean add(Integer i){
			//if(this.contains(i)) throw new RuntimeException("error");
			return (super.add(i));
		}
		
		@Override
		public int hashCode(){
			int hash = this.stream() .reduce(0, Integer::sum);
			return hash;
		}
		
		void roundBreaks(){
	   		for(int i =0; i<size(); i++){
	   			set(i,  TranscriptUtils.round(get(i), round));
	   		}
	   	}
		/*@Override
public boolean equals(Object o){
			CigarHash cc = (CigarHash) o;
		 return this.equals(cc.breaks);	
		}
	}*/
	}