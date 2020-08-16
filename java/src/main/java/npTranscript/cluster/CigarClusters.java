package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

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
	 this.refseq =null;
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

	CigarClusters( Sequence refSeq, int num_sources){
	//	this.thresh = thresh;
		this.refseq = refSeq;
		this.seqlen = refSeq.length();
		
		this.num_sources = num_sources;
	}

	public void update(Annotation annot){
		this.annot = annot;
	}
	
	Map<CigarHash, CigarCluster> l = new HashMap<CigarHash, CigarCluster>();

	
	
	public void clear() {
		l.clear();
		
	}
	int rem_count =0;
	public void matchCluster(CigarCluster c1,  int source_index, int num_sources, int chrom_index,String[] clusterIDs) throws NumberFormatException{
		String clusterID;
		CigarHash2 subID ;
		if(!l.containsKey(c1.breaks_hash)){
			CigarCluster newc = new CigarCluster("ID"+chrom_index+"."+(l.keySet().size()+rem_count), num_sources, c1, source_index);
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

	
	Annotation annot;
	final Sequence  refseq;
	final int seqlen;
	final int num_sources;
	public void process(CigarCluster cc,   Outputs o,String chrom,int chrom_index){
		cc.addZeros(seqlen); 
		String id = cc.id();
		int totalDepth = cc.readCountSum();
		int[] rcount=cc.readCount;
		
		if(o.clusterW!=null && totalDepth>IdentityProfile1.writeCoverageDepthThresh){
			int[][] matr =cc.getClusterDepth(num_sources, this.refseq);
			o.writeIntMatrix(id, matr);
		}
		if(o.altT!=null){
			if(IdentityProfile1.writeIsoformDepthThresh.length==1){
				if(totalDepth>=IdentityProfile1.writeIsoformDepthThresh[0]){
					o.writeString(id, cc.all_breaks, this.num_sources);
				}
			}else{
				boolean write=false;
				for(int j=0; j<rcount.length; j++){
					if(rcount[j]>=IdentityProfile1.writeIsoformDepthThresh[j]){
						write=true;
					}
				}
				if(write){
					o.writeString(id, cc.all_breaks, this.num_sources);
				}
			}
		}
		
	}
	public void process1(CigarCluster cc,   Outputs o,String chrom,int chrom_index, boolean forward,SortedSet<String> geneNames ){
		String read_count = TranscriptUtils.getString(cc.readCount);
		boolean hasLeaderBreak = TranscriptUtils.coronavirus  ? (cc.breaks.size()>1 &&  annot.isLeader(cc.breaks.get(1)*CigarHash2.round)) : false;
		geneNames.clear();
		String geneNme = annot.getString(cc.span, geneNames);
		o.printTranscript(
			cc.id()+"\t"+chrom+"\t"+cc.start+"\t"+cc.end+"\t"+annot.getTypeNme(cc.start, cc.end, forward)+"\t"+
	
		cc.exonCount()+"\t"+cc.numBreaks()+"\t"+(hasLeaderBreak? 1: 0)+"\t"+cc.breaks_hash.secondKey+"\t"+geneNme+"\t"+
		geneNames.size()+"\t"+
		cc.totLen+"\t"+cc.readCountSum()+"\t"+read_count
				+"\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false)+"\t"+cc.getErrorRatioSt());
	}
	
	/** clears consensus up to certain start position.  This designed to keep memory foot print under control.  Assumes the bams are sorted */
	public int clearUpTo(int endThresh, 	Outputs o, String chrom, int chrom_index,SortedSet<String> geneNames){
		if(this.l.size()< IdentityProfile1.clearUpThreshold) return 0;
		List<CigarHash> torem = new ArrayList<CigarHash>();
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();){
			CigarCluster cc = it.next();
			if(cc.end<endThresh){
				this.process1(cc, o, chrom, chrom_index, cc.forward, geneNames);
				this.process(cc, o, chrom, chrom_index);
				torem.add(cc.breaks_hash);
			}
		}
		//System.err.println("removed "+torem.size() + " of "+l.size());
		rem_count+=torem.size();
		for(int i=0; i<torem.size(); i++){
			l.remove(torem.get(i));
		}
		return torem.size();
	}
	
	
	public void getConsensus(  
			Outputs o, String chrom, int chrom_index,SortedSet<String> geneNames
			) throws IOException{
		if(TranscriptUtils.writeAnnotP) this.annot.print(o.annotP);
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
			CigarCluster cc = it.next();
			this.process1(cc, o, chrom, chrom_index, cc.forward, geneNames);
			this.process(cc, o, chrom, chrom_index);
		}
		
	}

	

}