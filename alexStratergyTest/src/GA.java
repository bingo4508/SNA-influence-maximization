import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

public class GA {
	public enum SelectionModel{RW_SELECT, TOURNAMENT_SELECT};
	public static MyRand myRand = new MyRand();

	public GA (SelectionModel n_sModel, FitnessFunction f, int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
		double n_pm, int n_timeLimit, int n_maxFe, int k, int round, Test.LTVertex[] arrVertices, int saveChromosomeNum, boolean firstRound)
	{
		//never = true;
		startTime = (int)System.currentTimeMillis()/1000;
		pEndTime = startTime;
		this.k = k;
		this.round = round;
		this.firstRound = firstRound;
		this.saveChromosomeNum = saveChromosomeNum;
		this.arrVertices = arrVertices;
		
		init (n_sModel, f, n_ell, n_nInitial, n_selectionPressure, n_pc, n_pm, n_timeLimit, n_maxFe);
	}

	public void init (SelectionModel n_sModel, FitnessFunction f, int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
			double n_pm, int n_timeLimit, int n_maxFe)
	{
		int i;

		sModel = n_sModel;
	    ell = n_ell;
	    nInitial = n_nInitial;
	    nCurrent = nInitial;
	    selectionPressure = n_selectionPressure;
	    pc = n_pc;
	    pm = n_pm;
	    timeLimit = n_timeLimit;
	    maxFe = n_maxFe;
	    stFitness = new Statistics();
	    

	    population = new Chromosome[nInitial];
	    offspring = new Chromosome[nInitial];
	    selectionIndex = new int[nInitial];

	    for (i = 0; i < nInitial; i++) {
	    	population[i] = new Chromosome(f, k, round, arrVertices);
	    	offspring[i] = new Chromosome(f, k, round, arrVertices);
	    }

	    initializePopulation ();
	}

	public void initializePopulation ()
	{
		int array[] = new int[arrVertices.length];
		BufferedReader input = null;
		String line = null;
		String[] nodes;
		int i = 0;
		if (!firstRound)
		{
			try 
			{
				input = new BufferedReader(new FileReader("GA_preSelectedGene.txt"));
				int scn = 0;
				for (scn = 0; scn < saveChromosomeNum && line != null; scn++)
					for (int r = 0; r < round && line != null; r++)
					{
						line = input.readLine();
						if (line != null)
						{
							nodes = line.split(" ");
							for (String node:nodes)
								population[scn].gene.get(r).add(Integer.parseInt(node));
						}
					}
				i = scn;
			} catch (FileNotFoundException f)
			{
				i = 0;
			}
			catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		for (; i < nInitial; i++)
		{
			myRand.uniformArray(array, 0, array.length, 0, array.length-1);
			for (int r = 0; r < round; r++)
			{
				for (int j = 0; j < k; j++)
					population[i].gene.get(r).add(array[r*k+j]);
			}
		}
	}

	public void selection ()
	{
		switch (sModel)
		{
			case RW_SELECT:
				rwSelection ();
				break;
			case TOURNAMENT_SELECT:
	    		tournamentSelection ();
	    		break;
	    	default:
	    		break;
	    }
	}

	/** tournament selection without replacement */
	public void tournamentSelection ()
	{
		int i, j;

	    // Adjusting population size 
	    nNextGeneration = getNextPopulation ();

	    if (nNextGeneration != nCurrent) {
	        selectionIndex = new int[nNextGeneration];
	        offspring = new Chromosome[nNextGeneration];

	        for (i = 0; i < nNextGeneration; i++)
	            offspring[i].init (k);
	    }

	    int randArray[] = new int[selectionPressure * nNextGeneration];

	    int q = (selectionPressure * nNextGeneration) / nCurrent;
	    int r = (selectionPressure * nNextGeneration) % nCurrent;

	    for (i = 0; i < q; i++)
	        myRand.uniformArray (randArray, (i * nCurrent), nCurrent, 0, nCurrent - 1);

	    myRand.uniformArray (randArray, (q * nCurrent), r, 0, nCurrent - 1);

	    for (i = 0; i < nNextGeneration; i++) {

	        int winner = 0;
	        double winnerFitness = -Double.MAX_VALUE;

	        for (j = 0; j < selectionPressure; j++) {
	            int challenger = randArray[selectionPressure * i + j];
	            double challengerFitness = population[challenger].getFitness ();

	            if (challengerFitness > winnerFitness) {
	                winner = challenger;
	                winnerFitness = challengerFitness;
	            }

	        }
	        selectionIndex[i] = winner;
	    }
	}

	/** Roulette wheel selection */
	public void rwSelection ()
	{
		int i, j;

	    // Adjusting population size 
	    nNextGeneration = getNextPopulation ();

	    if (nNextGeneration != nCurrent) {
	        selectionIndex = new int[nNextGeneration];
	        offspring = new Chromosome[nNextGeneration];

	        for (i = 0; i < nNextGeneration; i++)
	            offspring[i].init (k);
	    }

	    double totalFitness = 0.0;
	    for (i = 0; i < nCurrent; i++) 
		totalFitness += population[i].getFitness();

	    for (i = 0; i < nNextGeneration; i++) {
		double pointer = totalFitness * myRand.uniform();
		int index = -1;
		double partialSum = 0.0;
		for (j = 0; j < nCurrent; j++) {
		    partialSum += population[j].getFitness();
	            if (partialSum >= pointer) {
	                index = j;
	                break;
	            }
		}
		if (index == -1) index = nCurrent - 1;

		selectionIndex[i] = index;
	    }
	}

	public void crossover ()
	{
		if ((nNextGeneration & 0x1) == 0) { 
	    	// nNextGeneration is even
	    	
	        for (int i = 0; i < nNextGeneration; i += 2)
	        	mySNA_XO(population[selectionIndex[i]], population[selectionIndex[i + 1]],
	                offspring[i], offspring[i + 1]);
	
	    }
	    else {
	        for (int i = 0; i < nNextGeneration - 1; i += 2) {
	            mySNA_XO(population[selectionIndex[i]], population[selectionIndex[i + 1]],
	            	offspring[i], offspring[i + 1]);
	        }
	        offspring[nNextGeneration - 1] =
	            population[selectionIndex[nNextGeneration - 1]];
	    }
	}
	
	private void mySNA_XO(final Chromosome p1, final Chromosome p2, Chromosome c1, Chromosome c2)
	{
		for (int r = 0; r < round; r++)
		{
			HashSet<Integer> seedsInBothGene = new HashSet<Integer>();
			seedsInBothGene.addAll(p1.gene.get(r));
			seedsInBothGene.retainAll(p2.gene.get(r));
			HashSet<Integer> seedsOnlyInOne = new HashSet<Integer>();
			seedsOnlyInOne.addAll(p1.gene.get(r));
			seedsOnlyInOne.removeAll(seedsInBothGene);
			Integer seedIndex1[] = seedsOnlyInOne.toArray(new Integer[0]);
			seedsOnlyInOne.clear();
			seedsOnlyInOne.addAll(p2.gene.get(r));
			seedsOnlyInOne.removeAll(seedsInBothGene);
			Integer seedIndex2[] = seedsOnlyInOne.toArray(new Integer[0]);
			int randArray1[] = new int [seedIndex1.length];
			int randArray2[] = new int [seedIndex2.length];
			myRand.uniformArray(randArray1, 0, seedIndex1.length, 0, seedIndex1.length-1);
			myRand.uniformArray(randArray2, 0, seedIndex2.length, 0, seedIndex2.length-1);
			c1.gene.get(r).clear();
			c2.gene.get(r).clear();
			for (int i = 0; i < seedIndex1.length; i++)
			{
				if (myRand.uniform() < pc)
				{
					c1.gene.get(r).add(seedIndex2[randArray2[i]]);
					c2.gene.get(r).add(seedIndex1[randArray1[i]]);
				}
				else
				{
					c1.gene.get(r).add(seedIndex1[randArray1[i]]);
					c2.gene.get(r).add(seedIndex2[randArray2[i]]);
				}
			}
			for (Integer ind:seedsInBothGene)
			{
				c1.gene.get(r).add(ind);
				c2.gene.get(r).add(ind);
			}
		}
	}

	public void mutation ()
	{
		simpleMutation ();
	    //mutationClock ();
	}
	
	public void simpleMutation()
	{
		int q;
		HashSet<Integer> allIndex = new HashSet<Integer>(arrVertices.length);
	    for (int i = 0; i < arrVertices.length; i++)
	    {
	    	allIndex.add(i);
	    }
	    
		for (Chromosome c:offspring)
		{
			if (myRand.uniform()<pm)
			{
				q = myRand.uniformInt(0, ell-1);
				
				HashSet<Integer> seedsIndex = new HashSet<Integer> (k);
				HashSet<Integer> notSeedsIndex = new HashSet<Integer> (arrVertices.length-k);
				for (int r = 0; r < round; r++)
					seedsIndex.addAll(c.gene.get(r));
				notSeedsIndex.addAll(allIndex);
				notSeedsIndex.removeAll(seedsIndex);
				
				if (seedsIndex.contains(q))
				{
					for (int r = 0; r < round; r++)
						if (c.gene.get(r).contains(q))
						{
							c.gene.get(r).remove(q);
							Integer arrNotSeedsIndex[] = notSeedsIndex.toArray(new Integer[0]);
							c.gene.get(r).add(arrNotSeedsIndex[myRand.uniformInt(0, arrNotSeedsIndex.length-1)]);
							break;
						}
				}
				else
				{
					int r = myRand.uniformInt(0, round-1);
					Integer arrSeedsIndex[] = c.gene.get(r).toArray(new Integer[0]);
					c.gene.get(r).add(q);
					c.gene.get(r).remove(arrSeedsIndex[myRand.uniformInt(0, arrSeedsIndex.length-1)]);
				}
			}
		}
	}
	
	/*public void mutationClock ()
	{
		if (pm <= 1e-6) return; // can't deal with too small pm

	    int pointer = (int) (Math.log(1-myRand.uniform()) / Math.log(1-pm) + 1);
	    
	    HashSet<Integer> allIndex = new HashSet<Integer>(arrVertices.length);
	    for (int i = 0; i < arrVertices.length; i++)
	    {
	    	allIndex.add(i);
	    }
	    while (pointer < nNextGeneration * ell) {

		int q = pointer / ell;
		int r = pointer % ell;

		HashSet<Integer> seedsIndex = new HashSet<Integer> (k);
		HashSet<Integer> notSeedsIndex = new HashSet<Integer> (arrVertices.length-k);
		seedsIndex.addAll(offspring[q].gene);
		notSeedsIndex.addAll(allIndex);
		notSeedsIndex.removeAll(seedsIndex);
		
		if (offspring[q].gene.contains(r))
		{
			offspring[q].gene.remove(r);
			Integer arrNotSeedsIndex[] = notSeedsIndex.toArray(new Integer[0]);
			offspring[q].gene.add(arrNotSeedsIndex[myRand.uniformInt(0, arrNotSeedsIndex.length-1)]);
		}
		else
		{
			offspring[q].gene.add(r);
			Integer arrSeedsIndex[] = seedsIndex.toArray(new Integer[0]);
			offspring[q].gene.remove(arrSeedsIndex[myRand.uniformInt(0, arrSeedsIndex.length-1)]);
		}
		// Compute next mutation clock
		pointer += (int) (Math.log(1-myRand.uniform()) / Math.log(1-pm) + 1);
	    };
	}*/

	public void replacePopulation ()
	{
		int i;

	    if (nNextGeneration != nCurrent) {
	        population = new Chromosome[nNextGeneration];
	    }

	    for (i = 0; i < nNextGeneration; i++)
	        population[i].copyFrom(offspring[i]);

	    nCurrent = nNextGeneration;
	}

	//public boolean never;
	public void showStatistics ()
	{
		System.out.printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n",
		        generation, stFitness.getMax (), stFitness.getMean (),
		        stFitness.getMin ());
	    //System.out.printf ("best chromosome:");
	    //population[bestIndex].printf ();
	    System.out.printf ("\n");
	    /*if (population[bestIndex].fitness > 880 && never)
	    {
	    	Test.output.print("Get 800, "+generation+", ");
	    	never = false;
	    }*/
	}
	
	public void oneRun (boolean output)
	{
		int i;

		selection ();
	    crossover ();
	    mutation ();
	    replacePopulation ();
	    
	    double max = -Double.MAX_VALUE;
	    stFitness.reset ();
	    for (i = 0; i < nCurrent; i++) {
	        double fitness = population[i].getFitness ();
	        if (fitness > max) {
	            max = fitness;
	            bestIndex = i;
	        }
	        stFitness.record (fitness);
	    }

	    if (output)
	        showStatistics ();

	    generation++;
	}
	
	public void oneRun()
	{
		oneRun(true);
	}
	
	public int doIt (boolean output)
	{
		generation = 0;
		
	    while (!shouldTerminate ()) {
	        oneRun (output);
	    }
	    /*if (never)
	    {
	    	Test.output.print("Not get 800, "+generation+", ");
	    }*/
	    return generation;
	}
	
	public int doIt ()
	{
		return doIt(true);
	}

	public boolean shouldTerminate ()
	{
		boolean termination = false;

	    // Reach maximal # of function evaluations
	    if (maxFe != -1) {
	        if (fe > maxFe)
	            termination = true;
	    }

	    // Reach timeLimit - 15s
	    if (timeLimit != -1) {
	    	int endTime = (int)System.currentTimeMillis()/1000; 
	        if (endTime-startTime > timeLimit-(endTime-pEndTime)-15)
	            termination = true;
	        pEndTime = endTime;
	    }

	    // Found a satisfactory solution
	    if (stFitness.getMax() >= population[0].getMaxFitness())
	        termination = true;

	    // The population loses diversity
	    if (stFitness.getMax()-1e-6 < stFitness.getMean())
			termination = true;

	    return termination;
	}
	
	public int getNextPopulation ()
	{
		return nCurrent;
	}
	
	public static void outputErrMsg (final String errMsg)
	{
	    System.out.printf("%s\n", errMsg);
	    return;
	}
	
	public Chromosome getBestChromosome ()
	{
		return population[bestIndex];
	}
        

	public Statistics stFitness;

    protected int ell;                 // chromosome length
    protected int k;
    protected int round;
    protected boolean firstRound;
    protected int saveChromosomeNum;
    protected int nInitial;            // initial population size
    protected int nCurrent;            // current population size
    protected int nNextGeneration;     // population size for the next generation
    protected int selectionPressure;

    protected SelectionModel sModel;
    protected double pc;               // prob that XO happened when pairwise XO is used, 
								//or prob that each bit is flipped to be 1 when uniform XO is used
    protected double pm;               // prob of Mutation
    protected Chromosome []population;
    protected Chromosome []offspring;
    protected int []selectionIndex;
    protected int timeLimit;
    protected int maxFe;
    protected int repeat;
    protected int fe;
    protected int generation;
    protected int bestIndex;
    protected int startTime;
    protected int pEndTime;
    
    
    Test.LTVertex[] arrVertices;
}
