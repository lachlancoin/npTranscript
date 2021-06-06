package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import japsa.seq.Sequence;
import npTranscript.cluster.CigarCluster.Count;


/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
//	final double thresh;
	final int[] totalCounts;
	//final String chr;
	public CigarClusters(int num_sources, List<Sequence> genomes1){
	 this.num_sources = num_sources;
	// this.seqlen =0;
	// this.refseq =null;
	// this.annot = null;
	 this.totalCounts = new int[num_sources];
	 if(genomes1!=null){
		 this.genomes = new HashMap<String, Sequence>();
		 for(int i=0; i<genomes1.size(); i++){
			 this.genomes.put(genomes1.get(i).getName(), genomes1.get(i));
		 }
	 }else{
		 genomes = null;
	 }
	// chr = null;
	}
	
	private static String getKey(String[] str, int[] inds) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<inds.length; i++) {
			sb.append(str[inds[i]]);
			sb.append(",");
		}
		return sb.toString();
	}

	

	/*public void update(Annotation annot){
		this.annot = annot;
	}*/
	
	//SortedMap<CigarHash, CigarCluster> l = new TreeMap<CigarHash, CigarCluster>();
	Map<CigarHash, CigarCluster> l = new ConcurrentHashMap<CigarHash, CigarCluster>();
	
	
	/*public void clear() {
		l.clear();
		
	}*/
//	int rem_count =0;
	static String zero ="0";

	int currID=0;
	
	public synchronized void matchCluster(CigarCluster c1,  int source_index, int num_sources, String chr,String[] clusterIDs, 
			String strand, String readId) throws NumberFormatException{
		String clusterID;
		CigarHash br = c1.breaks_hash;
		CigarHash2 subID ;
		
		CigarCluster clust = l.get(br);
		this.totalCounts[source_index]++;
		if(clust==null){
			String id =  "ID"+"."+currID;
			currID++;
			String subID1 = 
					Outputs1.firstIsTranscriptome && source_index==0 ? readId:	id+".t0";
				
			CigarCluster newc = new CigarCluster(chr, id, num_sources, c1, source_index, strand, subID1);
			clusterID = newc.id();
			
			l.put(newc.breaks_hash, newc);
			subID = newc.breaks;
			clusterIDs[1] = subID1;
		}else{
			
			subID = clust.merge(c1, num_sources, source_index, clusterIDs, readId);
			clusterID = clust.id();
			
		}
	clusterIDs[0] = clusterID;
	//clusterIDs[1] = subID.rescale().toString();
		//return clusterID;
	}
	
	//final Annotation annot1;
	//final Sequence  refseq;
	//final int seqlen;
	final int num_sources;
	
	public static String getString(Object[] l){
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<l.length; i++){
			if(i>0) sb.append(",");
			sb.append(l[i].toString());
		}
		return sb.toString();
	}
	
	/*public void infoString(CigarCluster cc, HDFObjT obj, SortedSet<String> geneNames, String chrom){
		//String chrom = seq.getName();
		//String chrom = "";
		boolean forward = cc.forward;
		boolean hasLeaderBreak = TranscriptUtils.coronavirus  ? (cc.breaks.size()>1 &&  annot.isLeader(cc.breaks.get(1)*CigarHash2.round)) : false;
		geneNames.clear();
		int type_ind = annot.getTypeInd(cc.start, cc.end, forward);
		String type_nme = annot.nmes[type_ind];
		String geneNme = annot.getString(cc.span, geneNames);
		//obj[0] = cc.id(); 
		obj.chrom=chrom; obj.start = cc.start; obj.end = cc.end; obj.type_nme = type_nme;
		obj.exon_count = cc.exonCount(); obj.span = geneNames.size(); obj.genes = geneNme;
			
		 //obj[4] =type_nme;
		//obj[5] = cc.exonCount(); obj[6] = cc.numIsoforms()+""; obj[7] = cc.breaks_hash.secondKey; obj[8] = geneNme;
		//obj[9] = geneNames.size()+""; obj[10] = cc.totLen+""; obj[11] = cc.readCountSum()+"";
		geneNames.clear();
		//return getString(obj);
	}*/
	
	final Map<String, Sequence> genomes;
	
	
	//CigarHash fromKey = new CigarHash("",0);
	/** clears consensus up to certain start position.  This designed to keep memory foot print under control.  Assumes the bams are sorted */
	public synchronized int clearUpTo(int endThresh, 	Outputs o){
		//fromKey.end = endThresh;
	//	SortedMap<CigarHash, CigarCluster> tm = ((SortedMap)l).headMap(fromKey);
	//	if(tm.size()==0) return 0;
		int rem_count =0;
	//	List<CigarHash> torem = new ArrayList<CigarHash>();
			//	SortedSet<String> geneNames = new TreeSet<String>();
			//	l.tailMap(fromKey)
			
				for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();){
					CigarCluster cc = it.next();
					if(cc.endPos<endThresh){
						//int totalDepth = cc.readCountSum();
						//cc.process1( o,genomes==null ? null : genomes.get(cc.chrom));//, chrom, chrom_index);//,  geneNames);
					//	if(Outputs.writeIsoforms) o.writeIsoforms(cc, this,  totalDepth);
						o.writeDepthH5(cc.breaks_hash.secondKey, cc.map);
						rem_count++;
						l.remove(cc.breaks_hash);
						//l.put(cc.breaks_hash, null); //set as null.  If we remove it creates a potential problem with threading
						
						//torem.add(cc.breaks_hash);
					}
				}
		/*rem_count+=torem.size();
		for(int i=0; i<torem.size(); i++){
			l.remove(torem.get(i));
		}*/
		return rem_count;//torem.size();
	}
	
	
	
	public Iterator<CigarCluster> iterator(Comparator comp){
		if(comp==null) return l.values().iterator();
		List<CigarCluster> ss = new ArrayList<CigarCluster>(l.values());
		Collections.sort(ss, comp);
		return ss.iterator();
	}
	
	
	public void annotate(Annotation annot, Map<String, Sequence> genomes) {
		Iterator<CigarCluster> it = this.iterator(compPos);
		String currChrom = "";
		int seqlen  = 0;
		Set<String> done = new HashSet<String>();
		while(it.hasNext()){
			CigarCluster nxt = it.next();
			if(!nxt.chrom.equals(currChrom)){
				if(!done.add(currChrom)) throw new RuntimeException(" already had ");
				currChrom = nxt.chrom;
				if(genomes!=null){
					seqlen = genomes.get(currChrom).length();
				}
				annot.updateChrom(currChrom);
			}
			annot.annotate(nxt);
		}
		annot.close();
	}
	static  Comparator<CigarCluster> comp = new Comparator<CigarCluster>(){

		@Override
		public int compare(CigarCluster o1, CigarCluster o2) {
			// TODO Auto-generated method stub
		return -1*Integer.compare(o1.readCountSum(), o2.readCountSum());
		}
		
	};
	static  Comparator<CigarCluster> compPos = new Comparator<CigarCluster>(){

		@Override
		public int compare(CigarCluster o1, CigarCluster o2) {
			// TODO Auto-generated method stub
		return o1.breaks_hash.secondKey.compareTo(o2.breaks_hash.secondKey);
		}
		
	};
	
	public void getConsensus(  
			Outputs o
			) throws IOException{
		
		o.writeTotalString(this);
		Iterator<CigarCluster> it
			 = iterator(TranscriptUtils.coronavirus  ? comp : null);
			while(  it.hasNext()) {
					CigarCluster nxt = it.next();
				//	nxt.process1(o,genomes==null ? null : genomes.get(nxt.chrom));//, chrom, chrom_index, geneNames);
					int  totalDepth = nxt.readCountSum();
				//	if(Outputs.writeIsoforms) o.writeIsoforms(nxt, this, totalDepth);
					o.writeDepthH5(nxt.breaks_hash.secondKey, nxt.map);
			}
			
	}

	

}