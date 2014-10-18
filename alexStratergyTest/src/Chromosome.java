import java.util.ArrayList;
import java.util.HashSet;


public class Chromosome {
	
	public Chromosome (FitnessFunction f, int k, int round, Test.LTVertex arrVertices[])
	{
		this.arrVertices = arrVertices;
		this.round = round;
		gene = new ArrayList<HashSet<Integer>>(round);
	    init (k);
	    this.f = f;
	}
	
	public Chromosome copyFrom(Chromosome c)
	{
	    evaluated = c.evaluated;
	    fitness = c.fitness;
	    arrVertices = c.arrVertices;

    	for (int i = 0; i < 10; i++)
    	{
    		gene.get(i).clear();
    		gene.get(i).addAll(c.gene.get(i));
    	}

	    return this;
	}

	public void init (int k)
	{
		this.k = k;
		for (int r = 0; r < round; r++)
		{
			gene.add(new HashSet<Integer>(k));
		}
	    evaluated = false;
	}

	public double getFitness ()
	{
		if (evaluated)
	        return fitness;
	    else
	        return (fitness = evaluate ());
	}

	/** real evaluator */
	public double evaluate ()
	{
		evaluated = true;
	    return f.evaluate(this);
	}

	public boolean isEvaluated ()
	{
		return evaluated;
	}

	public void printf ()
	{
		for (int r = 0; r < round; r++)
		{
			System.out.printf("round"+(r+1)+": ");
		    for (Integer v:gene.get(r))
		        System.out.printf ("%d ", arrVertices[v].vertexIndex);
			System.out.println();
		}
	}


	public double getMaxFitness ()
	{
		// For OneMax
	    //return ((double)length-1e-6);
		return f.maxFitness();
	}

	public ArrayList<HashSet<Integer>> gene;
	protected int k;
	protected int round;
	protected double fitness;
	protected boolean evaluated;
	protected FitnessFunction f;
	protected Test.LTVertex[] arrVertices; 
}
