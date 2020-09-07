package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import japsa.seq.Sequence;


/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
//	final double thresh;
	
	public CigarClusters(int num_sources){
	 this.num_sources = num_sources;
	 this.seqlen =0;
	// this.refseq =null;
	 this.annot = null;
	}
	
	private static String getKey(String[] str, int[] inds) {
		StringBuffer sb = new StringBuffer();
		for(int i=0; i<inds.length; i++) {
			sb.append(str[inds[i]]);
			sb.append(",");
		}
		return sb.toString();
	}

	CigarClusters( Sequence refSeq, int num_sources, Annotation annot){
	//	this.thresh = thresh;
	//	this.refseq = refSeq;
		this.seqlen = refSeq.length();
		this.annot = annot;
		this.num_sources = num_sources;
	}

	/*public void update(Annotation annot){
		this.annot = annot;
	}*/
	
	SortedMap<CigarHash, CigarCluster> l = new TreeMap<CigarHash, CigarCluster>();

	
	
	/*public void clear() {
		l.clear();
		
	}*/
	int rem_count =0;
	
	public synchronized void matchCluster(CigarCluster c1,  int source_index, int num_sources, int chrom_index,String[] clusterIDs, char strand) throws NumberFormatException{
		String clusterID;
		CigarHash2 subID ;
		if(!l.containsKey(c1.breaks_hash)){
			CigarCluster newc = new CigarCluster("ID"+chrom_index+"."+(l.keySet().size()+rem_count), num_sources, c1, source_index, strand);
			clusterID = newc.id();
			l.put(newc.breaks_hash, newc);
			subID = newc.breaks;
		}else{
			CigarCluster clust = l.get(c1.breaks_hash);
			subID = clust.merge(c1, num_sources, source_index);
			clusterID = clust.id();
		}
	clusterIDs[0] = clusterID;
	
	clusterIDs[1] = subID.toString(CigarHash2.round, 1, subID.size()-1).hashCode()+"";
		//return clusterID;
	}

	
	final Annotation annot;
	//final Sequence  refseq;
	final int seqlen;
	final int num_sources;
	public void process1(CigarCluster cc,   Outputs o,Sequence seq, int chrom_index, SortedSet<String> geneNames ){
		String read_count = TranscriptUtils.getString(cc.readCount);
		String chrom = seq.getName();
		boolean forward = cc.forward;
		boolean hasLeaderBreak = TranscriptUtils.coronavirus  ? (cc.breaks.size()>1 &&  annot.isLeader(cc.breaks.get(1)*CigarHash2.round)) : false;
		geneNames.clear();
		int type_ind = annot.getTypeInd(cc.start, cc.end, forward);
		String type_nme = annot.nmes[type_ind];
		String geneNme = annot.getString(cc.span, geneNames);
	if(Outputs.writeGFF){
		cc.writeGFF(o.gffW, o.refOut[type_ind], chrom,  Outputs.isoThresh, type_nme, seq);
	}
		o.printTranscript(
			cc.id()+"\t"+chrom+"\t"+cc.start+"\t"+cc.end+"\t"+type_nme+"\t"+
	
		cc.exonCount()+"\t"+cc.numIsoforms()+"\t"+(hasLeaderBreak? 1: 0)+"\t"+cc.breaks_hash.secondKey+"\t"+geneNme+"\t"+
		geneNames.size()+"\t"+
		cc.totLen+"\t"+cc.readCountSum()+"\t"+read_count,"");
//				CigarCluster.recordDepthByPosition ?  cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false): "");
	}
	CigarHash fromKey = new CigarHash("",0);
	/** clears consensus up to certain start position.  This designed to keep memory foot print under control.  Assumes the bams are sorted */
	public synchronized int clearUpTo(int endThresh, 	Outputs o, Sequence chrom, int chrom_index){
		fromKey.end = endThresh;
		SortedMap<CigarHash, CigarCluster> tm = l.headMap(fromKey);
		if(tm.size()==0) return 0;
		List<CigarHash> torem = new ArrayList<CigarHash>();
				SortedSet<String> geneNames = new TreeSet<String>();
			//	l.tailMap(fromKey)
			
				for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();){
					CigarCluster cc = it.next();
					if(cc.end<endThresh){
						int totalDepth = cc.readCountSum();
						process1(cc, o, chrom, chrom_index,  geneNames);
						o.writeIsoforms(cc, this, chrom, chrom_index, totalDepth);
						o.writeDepthH5(cc, this, chrom,chrom_index, totalDepth);
						
						
						torem.add(cc.breaks_hash);
					}
				}
		rem_count+=torem.size();
		for(int i=0; i<torem.size(); i++){
			l.remove(torem.get(i));
		}
		return torem.size();
	}
	
	
	public void getConsensus(  
			Outputs o, Sequence chrom, int chrom_index,SortedSet<String> geneNames
			) throws IOException{
		if(TranscriptUtils.writeAnnotP) {
			Runnable run = new Runnable(){
				public void run(){
					annot.print(o.annotP);
				}
			};
			Outputs.h5writer.execute(run);
		}
		
			for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
					CigarCluster nxt = it.next();
					process1(nxt, o, chrom, chrom_index, geneNames);
					int  totalDepth = nxt.readCountSum();
					o.writeIsoforms(nxt, this, chrom, chrom_index, totalDepth);
					o.writeDepthH5(nxt, this, chrom,chrom_index, totalDepth);
			}
			
	}

	

}