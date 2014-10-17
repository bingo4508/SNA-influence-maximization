import java.util.HashSet;


public class Chromosome {
	
	public Chromosome (FitnessFunction f, int k, Test.LTVertex arrVertices[])
	{
		this.arrVertices = arrVertices;
		gene = null;
	    init (k);
	    this.f = f;
	}
	
	public Chromosome copyFrom(Chromosome c)
	{
	    evaluated = c.evaluated;
	    fitness = c.fitness;
	    arrVertices = c.arrVertices;
	    
	    gene.clear();
	    gene.addAll(c.gene);

	    return this;
	}

	public void init (int k)
	{
		this.k = k;
	    gene = new HashSet<Integer>(k);
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
	    for (Integer v:gene)
	        System.out.printf ("%d ", arrVertices[v].vertexIndex);
	}


	public double getMaxFitness ()
	{
		// For OneMax
	    //return ((double)length-1e-6);
		return f.maxFitness();
	}

	public HashSet<Integer> gene;
	protected int k;
	protected double fitness;
	protected boolean evaluated;
	protected FitnessFunction f;
	protected Test.LTVertex[] arrVertices; 
}
