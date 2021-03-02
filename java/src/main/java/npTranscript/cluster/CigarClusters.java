package npTranscript.cluster;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

import japsa.seq.Sequence;


/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
//	final double thresh;
	final int[] totalCounts;
	final String chr;
	public CigarClusters(int num_sources){
	 this.num_sources = num_sources;
	 this.seqlen =0;
	// this.refseq =null;
	 this.annot = null;
	 this.totalCounts = new int[num_sources];
	 chr = null;
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
		this.totalCounts = new int[num_sources];
		chr = refSeq.getName();
	}

	/*public void update(Annotation annot){
		this.annot = annot;
	}*/
	
	//SortedMap<CigarHash, CigarCluster> l = new TreeMap<CigarHash, CigarCluster>();
	Map<CigarHash, CigarCluster> l = new ConcurrentHashMap<CigarHash, CigarCluster>();
	
	
	/*public void clear() {
		l.clear();
		
	}*/
	int rem_count =0;
	
	public synchronized void matchCluster(CigarCluster c1,  int source_index, int num_sources, int chrom_index,String[] clusterIDs, char strand) throws NumberFormatException{
		String clusterID;
		CigarHash br = c1.breaks_hash;
		CigarHash2 subID ;
		CigarCluster clust = l.get(br);
		this.totalCounts[source_index]++;
		if(clust==null){
			
			CigarCluster newc = new CigarCluster("ID"+chrom_index+"."+(l.keySet().size()+rem_count), num_sources, c1, source_index, strand);
			clusterID = newc.id();
			l.put(newc.breaks_hash, newc);
			subID = newc.breaks;
			clusterIDs[1] = 0+"";
		}else{
			
			subID = clust.merge(c1, num_sources, source_index, clusterIDs);
			clusterID = clust.id();
			
		}
	clusterIDs[0] = clusterID;
	//clusterIDs[1] = subID.rescale().toString();
		//return clusterID;
	}
	
	final Annotation annot;
	//final Sequence  refseq;
	final int seqlen;
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
	
	public void process1(CigarCluster cc,   Outputs o,Sequence seq, int chrom_index, SortedSet<String> geneNames ){
		String read_count = TranscriptUtils.getString(cc.readCount);
		String chrom = seq.getName();
		boolean forward = cc.forward;
		boolean hasLeaderBreak = TranscriptUtils.coronavirus  ? (cc.breaks.size()>1 &&  annot.isLeader(cc.breaks.get(1)*CigarHash2.round)) : false;
		geneNames.clear();
		int type_ind = annot.getTypeInd(cc.startPos, cc.endPos, forward);
		String type_nme = annot.nmes[type_ind];
		String geneNme = annot.getString(cc.span, geneNames);
	if(Outputs.writeGFF){
		cc.writeGFF(o.gffW, o.refOut[type_ind], o.bedW,chrom,  Outputs.isoThresh, type_nme, seq);
	}
	String depth_str = "\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false);
		
		o.printTranscriptAlt(cc);
		o.printTranscript(
			cc.id()+"\t"+chrom+"\t"+cc.startPos+"\t"+cc.endPos+"\t"+type_nme+"\t"+
	
		cc.exonCount()+"\t"+cc.numIsoforms()+"\t"+(hasLeaderBreak? 1: 0)+"\t"+cc.breaks_hash.secondKey+"\t"+geneNme+"\t"+
		geneNames.size()+"\t"+
		cc.totLen+"\t"+cc.readCountSum()+"\t"+read_count,depth_str);
		List<Integer> br = cc.getBreaks();
	
		StringBuffer starts = new StringBuffer();
		StringBuffer ends = new StringBuffer();
		StringBuffer chroms = new StringBuffer();
		StringBuffer strand = new StringBuffer();
			int len =0;
			if(br!=null){
			for(int i=0; i<br.size(); i+=2){
				len += br.get(i+1) - br.get(i);
				String ext = i<br.size()-2 ? ";" : "";
				starts.append(br.get(i)+ext); ends.append(br.get(i+1)+ext);
				strand.append("+"+ext); chroms.append(chrom+ext);
			}
			}else{
				starts.append(cc.startPos); ends.append(cc.endPos); chroms.append(cc.endPos); strand.append("+");
			}
			
		
		o.printFC(cc.breaks_hash.secondKey+"\t"+chroms.toString()+"\t"+starts.toString()+"\t"+ends.toString()+"\t"+"+;+\t"+len+"\t"+read_count);
	//	Geneid  Chr     Start   End     Strand  Length  /DataOnline/Data/Projects/corona_invitro/host_analysis/direct_cDNA/vero/vero_24hpi/merged/genome/infected/mo
	//	ID0.0   MT007544.1;MT007544.1   14;27385        66;29860        +;+     2529    0       0       0       0       0       0
		//I
//				CigarCluster.recordDepthByPosition ?  cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false): "");
	}
	//CigarHash fromKey = new CigarHash("",0);
	/** clears consensus up to certain start position.  This designed to keep memory foot print under control.  Assumes the bams are sorted */
	public synchronized int clearUpTo(int endThresh, 	Outputs o, Sequence chrom, int chrom_index){
		//fromKey.end = endThresh;
	//	SortedMap<CigarHash, CigarCluster> tm = ((SortedMap)l).headMap(fromKey);
	//	if(tm.size()==0) return 0;
		int rem_count =0;
	//	List<CigarHash> torem = new ArrayList<CigarHash>();
				SortedSet<String> geneNames = new TreeSet<String>();
			//	l.tailMap(fromKey)
			
				for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();){
					CigarCluster cc = it.next();
					if(cc.endPos<endThresh){
						int totalDepth = cc.readCountSum();
						process1(cc, o, chrom, chrom_index,  geneNames);
						if(Outputs.writeIsoforms) o.writeIsoforms(cc, this, chrom, chrom_index, totalDepth);
						o.writeDepthH5(cc, this, chrom,chrom_index, totalDepth);
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
	
	
	public void getConsensus(  
			Outputs o, Sequence chrom, int chrom_index,SortedSet<String> geneNames
			) throws IOException{
		if(TranscriptUtils.writeAnnotP) {
					annot.print(o.annotP);
		}
		o.writeTotalString(this);
		Iterator<CigarCluster> it;
		if(TranscriptUtils.coronavirus ){
		Comparator<CigarCluster> comp = new Comparator<CigarCluster>(){

			@Override
			public int compare(CigarCluster o1, CigarCluster o2) {
				// TODO Auto-generated method stub
			return -1*Integer.compare(o1.readCountSum(), o2.readCountSum());
			}
			
		};
		List<CigarCluster> ss = new ArrayList<CigarCluster>(l.values());
		Collections.sort(ss, comp);
			it = ss.iterator();
		}else{
			it = l.values().iterator();
		}
			while(  it.hasNext()) {
					CigarCluster nxt = it.next();
					process1(nxt, o, chrom, chrom_index, geneNames);
					int  totalDepth = nxt.readCountSum();
					if(Outputs.writeIsoforms) o.writeIsoforms(nxt, this, chrom, chrom_index, totalDepth);
					o.writeDepthH5(nxt, this, chrom,chrom_index, totalDepth);
			}
			
	}

	

}