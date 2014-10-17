import java.nio.ByteBuffer;
import java.util.Random;


public class MyRand {
	private MT19937 mtRand;
	public MyRand ()
	{
		Random rand = new Random();
		int N = 624;
		byte init_key[] = new byte[N*Long.BYTES];
		long temp;

	    for (int i = 0; i < N; i++)
	    {
	    	temp = (long) rand.nextLong() * rand.nextInt();
	    	for (int j = 0; j < Long.BYTES; j++)
	    		init_key[4*i+j] = ByteBuffer.allocate(8).putLong(temp).array()[j];
	    }

	    mtRand = new MT19937(init_key);
	}
	
	public boolean flip (double prob)
	{
		return (uniform () < prob);
	}
	
	public boolean flip ()
	{
		return flip(0.5);
	}

	/** From [a,b) */
	public double uniform (double a, double b)
	{
		return uniform () * (b - a) + a;
	}

	/** From [0,1) */
	public double uniform ()
	{
		return mtRand.nextDouble();
	}

	/** Int From [a,b] */
	public int uniformInt (int a, int b)
	{
		return (a + (int) (uniform () * (b - a + 1)));
	}

	/** Generate a random array of size num, from [a,b] */
	public void uniformArray (int []array, int startIndex, int num, int a, int b)
	{
		int []base = new int[b - a + 1];
	    int i;
	    int r;

	    for (i = 0; i < b - a + 1; i++)
	        base[i] = a + i;

	    for (i = 0; i < num; i++) {
	        r = uniformInt (0, b - a - i);
	        array[startIndex+i] = base[r];
	        base[r] = base[b - a - i];
	    }
	}
}
