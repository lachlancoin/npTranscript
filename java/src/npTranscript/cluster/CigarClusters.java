package npTranscript.cluster;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import npTranscript.cluster.IdentityProfile1.Outputs;

/**
 * @author Lachlan Coin
 *
 */


public class CigarClusters {
	
	
	
//	final double thresh;
	
	CigarClusters( Sequence refSeq, Annotation annot, int num_sources){
	//	this.thresh = thresh;
		this.refseq = refSeq;
		this.seqlen = refSeq.length();
		this.annot = annot;
		this.num_sources = num_sources;
	}
	/*public DistanceMatrix getDistanceMatrix( PrintWriter pw){
		CigarCluster[] l1 = this.l.values().toArray(new CigarCluster[0]);
		int len  = l1.length;
		double[][] res = new double[len][];
		String[] labels = new String[len];
		for(int i=0; i<len; i++) {
			
			res[i] = new double[len];
			CigarCluster cc = l1[i];
			labels[i] = cc.id;
			res[i][i] =0;
			for(int j=0; j<i; j++) {
				CigarCluster cc_j = l1[j];
				double dist = 1-cc.similarity(cc_j);
				res[i][j] = dist;
				res[j][i] = dist;
			}
			
		}
		for(int i=0; i<len; i++) {
			pw.print(labels[i]+","+l1[i].index+",");
			pw.println(getString(res[i]));
		}
	
		IdGroup group = new  SimpleIdGroup(labels);
		DistanceMatrix dm = new DistanceMatrix(res, group);
		return dm;
	}*/
	Map<CigarHash, CigarCluster> l = new HashMap<CigarHash, CigarCluster>();


	
	public String matchCluster(CigarCluster c1,  int source_index, int num_sources, String chrom) throws NumberFormatException{
	
		String clusterID = chrom;
	
		
		if(!l.containsKey(c1.breaks)){
			CigarCluster newc = new CigarCluster(l.keySet().size(), num_sources, c1, source_index);
			
			clusterID = newc.id(chrom);
			l.put(newc.breaks, newc);
			//System.err.println("new cluster " +" "+newc.id+" "+index);

		}else{
			CigarCluster clust = l.get(c1.breaks);
			clust.merge(c1);
			clusterID = clust.id(chrom);
		}
		return clusterID;
	}

	/*private void writeSeq(SequenceOutputStream seqFasta,Annotation annot,  PrintWriter exonP, int[][] exons, CigarCluster cc, String refseq, int[] firstlast){
		StringBuffer descline; String annotline;
		String id = cc.id+"";
		int read_count = cc.readCountSum;
		if(cc.breakSt>0){
			String leftseq = refseq.subSequence(cc.start, cc.breakSt).toString();
			String rightseq = refseq.subSequence(cc.breakEnd, cc.end).toString();
		}else{
			String rightseq = refseq.subSequence(cc.start, cc.end).toString();
		}
		for(int j=0; j<exons.length; j++) {
			int start = exons[j][0];
			int end = exons[j][1];
			exonP.println(id+"\t"+start+"\t"+end+"\t"+read_count);

			
			annotline.append(annot.calcORFOverlap(start, end, first_last, transcript_len));

			int len = end-start+1;
			descline.append(";");
			descline.append(start); descline.append("-"); descline.append(end); descline.append(","); descline.append(len);
			subseq.append();
			
			transcript_len += len;
			//seqline.append(subseq.toString());
		}
		Sequence subseq1 = new Sequence(refseq.alphabet(),subseq.toString().toCharArray(), id);
	
		descline.append(" "); descline.append(annotline);
		subseq1.setDesc(descline.toString());
		subseq1.writeFasta(seqFasta); 
	}*/
	//static int tolerance = 10;
	final Annotation annot;
	final Sequence  refseq;
	final int seqlen;
	final int num_sources;
	public String process(CigarCluster cc,   Outputs o,
			Map<String, List<CigarHash>> geneToHash,String chrom){
		cc.addZeros(seqlen); 
		
		 int[][] matr =cc.getClusterDepth(num_sources);
		String id = cc.id(chrom);
		o.writeIntMatrix(id, matr);
		
		String read_count = TranscriptUtils.getString(cc.readCount);
		int startPos, endPos;
		if(!IdentityProfile1.annotByBreakPosition ){
			startPos = cc.start;
			endPos =cc.end;
			

		}else{
			startPos = cc.breakSt;
			endPos = cc.breakEnd;
		}
		String downstream = annot==null ? null : annot.nextDownstream(endPos);
		String upstream = annot==null ? null : annot.nextUpstream(startPos);
		if(upstream==null) upstream =  TranscriptUtils.round(startPos, CigarHash.round)+"";
		if(downstream==null) downstream =  TranscriptUtils.round(endPos, CigarHash.round)+"";
		o.printTranscript(
			id+"\t"+chrom+"\t"+cc.start+"\t"+cc.end+"\t"+cc.getTypeNme(seqlen)+"\t"+
		cc.breaks.toString()+"\t"+cc.breaks.hashCode()+"\t"+
		cc.breakSt+"\t"+cc.breakEnd+"\t"+
		upstream+"\t"+downstream+"\t"+
		cc.totLen+"\t"+cc.readCountSum+"\t"+read_count
				+"\t"+cc.getTotDepthSt(true)+"\t"+cc.getTotDepthSt(false)+"\t"+cc.getErrorRatioSt());
		return TranscriptUtils.round(cc.start, 100)+"."+downstream+"."+upstream+"."+TranscriptUtils.round(seqlen-cc.end, 100);
	}
	
	public void getConsensus(  
			Outputs o, String chrom
			) throws IOException{
		
        Map<String, List<CigarHash>> geneToHash = new HashMap<String, List<CigarHash>>();
		for(Iterator<CigarCluster> it = l.values().iterator(); it.hasNext();) {
			CigarCluster cc = it.next();
			String key = this.process(cc, o, geneToHash,  chrom);
			if(geneToHash!=null){
				List<CigarHash> l;
				if(!geneToHash.containsKey(key)) geneToHash.put(key, l = new ArrayList<CigarHash>());
				else  l = geneToHash.get(key);
				l.add(cc.breaks);
			}
		}
		
		
		
	}

}