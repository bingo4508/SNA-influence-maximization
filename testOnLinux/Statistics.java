
public class Statistics {
	public Statistics () {
	    reset ();
	}

	public void reset () {
		precision = 1e-6;
		min = Double.MAX_VALUE-1;
		second_min = Double.MAX_VALUE-2;
		max = -Double.MAX_VALUE+1;
		second_max = -Double.MAX_VALUE+2;
		sum = 0.0;
		variance = 0.0;
		number = 0;
		status = true;
	}

	public void record (double value) {
	
	    if (status == false)
	        return;
	
	    number++;
	    sum += value;
	    variance += value * value;
	    if (min > value + precision) {
	        second_min = min;
	        min = value;
	    }
	    if (max < value - precision) {
	        second_max = max;
	        max = value;
	    }
	}

	/** get the number of samples */
	public long getNumber () {
	    return number;
	}

	/** get mean */
	public double getMean () {
	    return sum / number;
	}
	
	/** get variance */
	public double getVariance () {
	    double mean = getMean ();
	    return variance / number - mean * mean;
	}
	
	/** get standard deviation */
	public double getStdev () {
	    return Math.sqrt (getVariance ());
	}
	
	public double getMin () {
	    return min;
	}
	
	public double getMax () {
	    return max;
	}
	
	public double getSecondMax () {
	    return second_max;
	}
	
	public double getSecondMin () {
	    return second_min;
	}
	
	public void turnOn () {
	    status = true;
	}
	
	public void turnOff () {
	    status = false;
	}
        
	private double precision;
	private double min;
	private double max;
	private double second_min;
	private double second_max;
	private double sum;
	private double variance;
	private long number;
	private boolean status;
}
