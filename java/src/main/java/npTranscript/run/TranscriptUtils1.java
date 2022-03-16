package npTranscript.run;

import java.util.ArrayDeque;
import java.util.Comparator;
import java.util.TreeSet;


public class TranscriptUtils1 {

	
	public static double[] slidingWindowMeans(byte[] nums, int k) {
		// TODO Auto-generated method stub
		double[] means = new double[nums.length - k + 1];
		ArrayDeque<Byte> deque = new ArrayDeque<>();
		double sum = 0.0D;
		for(int i = 0; i <= nums.length; i++) {
			if(i >= k) {
				means[i - k] = sum / k;
				sum -= deque.pop();
			} 
			if(i < nums.length) {
				deque.offer(nums[i]);
				sum += nums[i];
			} 
		} 
		return means;
	}

	public static double[] slidingWindowPercentiles(byte[] nums, int k, double q) {
		// TODO Auto-generated method stub
		double[] percentiles = new double[nums.length - k + 1];
		SlidingWindowPercentile percentileWindow = new SlidingWindowPercentile(nums, q);
		for(int i = 0; i <= nums.length; i++) {
			if(i >= k) {
				percentiles[i - k] = percentileWindow.getPercentile();
				percentileWindow.remove(i - k);
			} 
			if(i < nums.length)
				percentileWindow.add(i); 
		} 
		return percentiles;
	}

	public static double percentile(byte[] nums, int startPos, int len, double q) {
		// TODO Auto-generated method stub
		if(startPos + len > nums.length)
			throw new RuntimeException("!!!"); 
		SlidingWindowPercentile percentileWindow = new SlidingWindowPercentile(nums, q);
		for(int i = 0; i < len; i++)
			percentileWindow.add(startPos + i); 
		return percentileWindow.getPercentile();
	}

	public static double percentile(byte[] nums, double q) {
		// TODO Auto-generated method stub
		return percentile(nums, 0, nums.length, q);
	}

	public static double median(byte[] nums, int startPos, int len) {
		if(nums.length==0) return Double.NaN;
		// TODO Auto-generated method stub
		return percentile(nums, startPos, len, 0.5D);
	}

	public static double median(byte[] nums, double q) {
		// TODO Auto-generated method stub
		return percentile(nums, 0, nums.length, 0.5D);
	}

	public static double[] slidingWindowMedians(byte[] nums, int k) {
		return slidingWindowPercentiles(nums, k, 0.5D);
	}

	private static class SlidingWindowPercentile {
		private TreeSet<Integer> lower;
		private TreeSet<Integer> upper;
		private byte[] elements;
		private double percentile;

		public SlidingWindowPercentile(byte[] nums, double p) {
			if(p < 0.0D || p > 1.0D)
				throw new RuntimeException("!!!"); 
			elements = nums;
			Comparator<Integer> c = new Comparator<Integer>() {
				public int compare(Integer a, Integer b) {
					return elements[a] == elements[b] ? a - b : elements[a] < elements[b] ? -1 : 1;
				}
			};
			lower = new TreeSet<>(c);
			upper = new TreeSet<>(c);
			percentile = p;
		}

		public void add(int i) {
			if(lower.size() * (1.0D - percentile) <= upper.size() * percentile) {
				upper.add(i);
				int m = upper.first();
				upper.remove(m);
				lower.add(m);
			} else {
				lower.add(i);
				int m = lower.last();
				lower.remove(m);
				upper.add(m);
			} 
		}

		public boolean remove(int i) {
			return lower.remove(i) || upper.remove(i);
		}

		public double getPercentile() {
			if(lower.isEmpty() && upper.isEmpty())
				return Double.NaN; 
			if(lower.isEmpty())
				return elements[upper.first()]; 
			if(upper.isEmpty())
				return elements[lower.last()]; 
			int n1 = lower.size();
			int n2 = upper.size();
			int n = n1 + n2 - 1;
			double p1 = (double)(n1 - 1) / n;
			double p2 = (double)n1 / n;
			return elements[lower.last()]+(elements[upper.first()]-elements[lower.last()])*(percentile-p1)/(p2-p1);
		}
	}
}
