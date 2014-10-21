import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Stack;

import edu.uci.ics.jung.graph.DirectedSparseMultigraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseMultigraph;
import edu.uci.ics.jung.graph.util.Pair;


public class Test {
	
	private static class DoubleRef
	{
		public double value;
		public DoubleRef(double v)
		{
			value = v;
		}
	}
	
	private static class VertexWithData
	{
		public double data;
		public LTVertex vertex;
		public VertexWithData(LTVertex v, double d)
		{
			vertex = v;
			data = d;
		}
	}
	
	private static class SortedVertexList extends ArrayList<VertexWithData>
	{
	    /**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public boolean add(VertexWithData mt)
	    {
	        int index = Collections.binarySearch(this, mt, new Comparator<VertexWithData>(){
	        	public int compare(VertexWithData o1, VertexWithData o2)
	        	{
	        		return (int)(o1.data - o2.data);
	        	}
	        });
	        if (index < 0) index = ~index;
	        super.add(index, mt);
	        return true;
	    }
	}
	
	private static class SNA_FitnessFunction implements FitnessFunction {
		int playerID;
		Graph<LTVertex, LTEdge> g;
		HashMap<Integer, LTVertex> vertices;
		HashMap<LTVertex, Pair<Double>> ltHoldMap;
		int round;
		protected HashSet<Test.LTVertex> player1NewlyActiveVertices;
	    protected HashSet<Test.LTVertex> player2NewlyActiveVertices; 
	    protected HashSet<Test.LTVertex> player1ActiveVertices;
	    protected HashSet<Test.LTVertex> player2ActiveVertices;
	    
		public SNA_FitnessFunction(int playerID, Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, HashMap<LTVertex, Pair<Double>> ltHoldMap, int round, 
				HashSet<Test.LTVertex> player1NewlyActiveVertices, HashSet<Test.LTVertex> player2NewlyActiveVertices,
				HashSet<Test.LTVertex> player1ActiveVertices, HashSet<Test.LTVertex> player2ActiveVertices)
		{
			this.playerID = playerID;
			this.g = g;
			this.vertices = vertices;
			this.ltHoldMap = ltHoldMap;
			this.round = round;
			this.player1NewlyActiveVertices = player1NewlyActiveVertices;
			this.player2NewlyActiveVertices = player2NewlyActiveVertices;
			this.player1ActiveVertices = player1ActiveVertices;
			this.player2ActiveVertices = player2ActiveVertices;
		}
		@Override
		public double evaluate(Chromosome c)
		{
			HashSet<LTVertex> seeds1 = new HashSet<LTVertex>();
			HashSet<LTVertex> seeds2 = new HashSet<LTVertex>();
			HashSet<LTVertex> activeVertices1 = new HashSet<LTVertex>();
			HashSet<LTVertex> activeVertices2 = new HashSet<LTVertex>();
			Pair<HashSet<LTVertex>> newerActiveVertices = new Pair<HashSet<LTVertex>>(new HashSet<LTVertex>(), new HashSet<LTVertex>());
			HashSet<Vertex> activatedVertices1= new HashSet<Vertex>();
			HashSet<Vertex> activatedVertices2= new HashSet<Vertex>();
			
			resetState(vertices, player1NewlyActiveVertices, player2NewlyActiveVertices, player1ActiveVertices, player2ActiveVertices);
			resetHold(g, vertices, ltHoldMap);
			for (int r = 0; r < round; r++)
			{
				seeds1.clear();
				seeds1.addAll(newerActiveVertices.getFirst());
				seeds2.clear();
				seeds2.addAll(newerActiveVertices.getSecond());
				if (r == 0)
				{
					seeds1.addAll(player1NewlyActiveVertices);
					seeds2.addAll(player2NewlyActiveVertices);
				}
				if (playerID == 1)
					for (Integer ind:c.gene.get(r))
						seeds1.add(c.arrVertices[ind]);
				else if (playerID == 2)
					for (Integer ind:c.gene.get(r))
						seeds2.add(c.arrVertices[ind]);
				activeVertices1.clear();
				activeVertices2.clear();
				newerActiveVertices = null;
				
				newerActiveVertices = runDiffusion(g, vertices, seeds1, seeds2, 0, activeVertices1, activeVertices2, Model.LTModel);
				activatedVertices1.addAll(activeVertices1);
				activatedVertices2.addAll(activeVertices2);
				activatedVertices1.addAll(newerActiveVertices.getFirst());
				activatedVertices2.addAll(newerActiveVertices.getSecond());
			}
			if (playerID == 1)
				return activatedVertices1.size()-activatedVertices2.size();
			else
				return activatedVertices2.size()-activatedVertices1.size();
		}
		@Override
		public double maxFitness()
		{
			return 1000000;	//maxFitness is unknown
		}
	}
	
	public static class Edge
	{
	}
	
	public static class Vertex
	{
		public static enum State{INACTIVE, PLAYER1_NEWLY_ACTIVE, PLAYER2_NEWLY_ACTIVE, PLAYER1_ACTIVE, PLAYER2_ACTIVE}; 
		public int vertexIndex;
		public State state;
		public Vertex(int vertexIndex, State state)
		{
			this.vertexIndex = vertexIndex;
			this.state = state;
		}
	}
	
	public static class LTEdge extends Edge
	{
		public double ltInfluence;
		public LTEdge(double ltInfluence)
		{
			this.ltInfluence = ltInfluence;
		}
	}
	
	public static class LTVertex extends Vertex
	{
		public double ltThreshold;
		public double ltHold1;
		public double ltHold2;
		public LTVertex(int vertexIndex, double ltThreshold, State state)
		{
			super(vertexIndex, state);
			this.ltThreshold = ltThreshold;
			this.ltHold1 = 0;
			this.ltHold2 = 0;
		}
	}
	
	public static enum Model{LTModel};
	
	private static class InputFileCollect
	{
		public String edgesFileName;
		public String nodesFileName;
		public String statusFileName;
		public Model model;
		public InputFileCollect(String edgesFileName, String nodesFileName, String statusFileName)
		{
			this.edgesFileName = edgesFileName;
			this.nodesFileName = nodesFileName;
			this.statusFileName = statusFileName;
			this.model = Model.LTModel;
		}
	}
	
	public static void main(String[] args)
	{
		//testStratergyForOnePartyGA();
		testStratergyForMultiParty(args);
		//testStratergyForOnePartyGreedy();
	}
	
	private static void testStratergyForMultiParty(String []args)
	{
		if (args.length != 7)
		{
			System.out.println("Format: Test $(PLAYER_ID) $(NODES_FILE) $(EDGES_FILE) $(STATUS_FILE) $(NODES_NUM_PER_ITER) $(SELECTED_NODES_FILE) $(TIME_LIMIT_IN_SEC)");
			return;
		}
		int playerID = Integer.parseInt(args[0]);
		InputFileCollect inputFile = new InputFileCollect(args[2], args[1], args[3]);
		int k = Integer.parseInt(args[4]);
		String outputFileName = args[5];
		int timeLimit = (int)Double.parseDouble(args[6]);
		
		Graph<LTVertex, LTEdge> g;
		HashMap<Integer, LTVertex> vertices;
		HashSet<LTVertex> seeds = null;
		
		vertices = new HashMap<Integer, LTVertex>();
		g = createGraph(inputFile, vertices);
		
		if (playerID == 1)
		{
			seeds = stratergyForPlayer1(g, vertices, k, inputFile, timeLimit);
		}
		else if (playerID == 2)
		{
			seeds = stratergyForPlayer2(g, vertices, k, inputFile, timeLimit);
		}
		
		
		PrintStream output = null;
		try {
			output = new PrintStream(new File(outputFileName));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (LTVertex s:seeds)
		{
			output.print(""+s.vertexIndex+" ");
		}
		output.close();
		return;
	}
	
	private static HashSet<LTVertex> stratergyForPlayer1(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, int k, InputFileCollect inputFile, int timeLimit)
	{
		BufferedReader input = null;
		String line[] = new String[4];
		String nodes[];
		HashSet<LTVertex> player1ActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player2ActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player1NewlyActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player2NewlyActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player1PreviousActiveVertices = new HashSet<LTVertex>();	//actually, previous seeds
		HashSet<LTVertex> player2PreviousActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player1PreviousActivatedVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player2PreviousActivatedVertices = new HashSet<LTVertex>();
		int round = 10;
		try 
		{
			input = new BufferedReader(new FileReader(inputFile.statusFileName));
			line[0] = input.readLine();
			line[1] = input.readLine();
			line[2] = input.readLine();
			line[3] = input.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		while(line[0] != null)
		{
			round--;
			nodes = line[0].split(" ");
			player1PreviousActiveVertices.clear();
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					if (!player2ActiveVertices.contains(v))
					{
						player1ActiveVertices.add(v);
						player1PreviousActiveVertices.add(v);
					}
				}
			}
			nodes = line[1].split(" ");
			player2PreviousActiveVertices.clear();
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					if (!player1ActiveVertices.contains(v))
					{
						player2ActiveVertices.add(v);
						player2PreviousActiveVertices.add(v);
					}
				}
			}
			nodes = line[2].split(" ");
			player1PreviousActivatedVertices.clear();
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					player1ActiveVertices.add(v);
					player1PreviousActivatedVertices.add(v);
				}
			}
			nodes = line[3].split(" ");
			player2PreviousActivatedVertices.clear();
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					player2ActiveVertices.add(v);
					player2PreviousActivatedVertices.add(v);
				}
			}
			try {
				line[0] = input.readLine();
				line[1] = input.readLine();
				line[2] = input.readLine();
				line[3] = input.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		HashMap<LTVertex, Pair<Double>> ltHoldMap = null;
		if (round != 10)
			ltHoldMap = createLTHoldMap(vertices);
		player1ActiveVertices.removeAll(player1PreviousActivatedVertices);
		player1ActiveVertices.removeAll(player1PreviousActiveVertices);
		player2ActiveVertices.removeAll(player2PreviousActivatedVertices);
		player2ActiveVertices.removeAll(player2PreviousActiveVertices);
		resetState(vertices, player1PreviousActiveVertices, player2PreviousActiveVertices, player1ActiveVertices, player2ActiveVertices);
		resetHold(g, vertices, ltHoldMap);
		updateLTHold(g, vertices, player1PreviousActiveVertices, player2PreviousActiveVertices);
		player1ActiveVertices.addAll(player1PreviousActivatedVertices);
		player1ActiveVertices.addAll(player1PreviousActiveVertices);
		player2ActiveVertices.addAll(player2PreviousActivatedVertices);
		player2ActiveVertices.addAll(player2PreviousActiveVertices);
		ltHoldMap = createLTHoldMap(vertices);
		try {
			input.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int N1 = 900, s1 = 2, maxG = 25;
		double pc1 = 0.5625, pm1 = 0.12;
		
		return GA(1, g, vertices, ltHoldMap, player1NewlyActiveVertices, player2NewlyActiveVertices, player1ActiveVertices, player2ActiveVertices, k, round, N1, s1, pc1, pm1, timeLimit).get(0);
	}
	
	private static HashSet<LTVertex> stratergyForPlayer2(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, int k, InputFileCollect inputFile, int timeLimit)
	{
		BufferedReader input = null;
		String line[] = new String[4];
		String nodes[];
		HashSet<LTVertex> player1ActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player2ActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player1NewlyActiveVertices = new HashSet<LTVertex>();
		HashSet<LTVertex> player2NewlyActiveVertices = new HashSet<LTVertex>();
		int round = 10;
		try 
		{
			input = new BufferedReader(new FileReader(inputFile.statusFileName));
			line[0] = input.readLine();
			line[1] = input.readLine();
			line[2] = input.readLine();
			line[3] = input.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		while(line[1] != null)
		{
			round--;
			nodes = line[0].split(" ");
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					if (!player2ActiveVertices.contains(v))
						player1ActiveVertices.add(v);
				}
			}
			nodes = line[1].split(" ");
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					if (!player1ActiveVertices.contains(v))
						player2ActiveVertices.add(v);
				}
			}
			nodes = line[2].split(" ");
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					player1ActiveVertices.add(v);
				}
			}
			nodes = line[3].split(" ");
			for (String node:nodes)
			{
				if (!node.equals(""))
				{
					LTVertex v = vertices.get(Integer.parseInt(node));
					player2ActiveVertices.add(v);
				}
			}

			try {
				line[0] = input.readLine();
				line[1] = input.readLine();
				line[2] = input.readLine();
				line[3] = input.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		nodes = line[0].split(" ");			//seeds player1 chose in this round
		for (String node:nodes)
			player1NewlyActiveVertices.add(vertices.get(Integer.parseInt(node)));
		HashMap<LTVertex, Pair<Double>> ltHoldMap = null;
		if (round != 10)
			ltHoldMap = createLTHoldMap(vertices);
		try {
			input.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int N1 = 900, s1 = 2, maxG = 25;
		double pc1 = 0.5625, pm1 = 0.12;
		
		HashSet<LTVertex> seeds = GA(2, g, vertices, ltHoldMap, player1NewlyActiveVertices, player2NewlyActiveVertices, player1ActiveVertices, player2ActiveVertices, k, round, N1, s1, pc1, pm1, timeLimit).get(0);
		resetState(vertices, player1NewlyActiveVertices, player2NewlyActiveVertices, player1ActiveVertices, player2ActiveVertices);
		resetHold(g, vertices, ltHoldMap);
		updateLTHold(g, vertices, player1NewlyActiveVertices, seeds);
		return seeds;
	}
	//public static PrintStream output = null;
	private static void testStratergyForOnePartyGA()
	{
		InputFileCollect[] inputFiles = new InputFileCollect[]{
				//new InputFileCollect("partB_hepth_lt_edges.txt", "partB_hepth_lt_nodes.txt", ""),
				new InputFileCollect("partB_egofb_lt_edges.txt", "partB_egofb_lt_nodes.txt", ""),
		};
		
		/*try {
			output = new PrintStream(new File("TestGA_WithDifferentPC_PM2.cvs"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
		Graph<LTVertex, LTEdge> g;
		HashMap<Integer, LTVertex> vertices;
		ArrayList<HashSet<LTVertex>> seedsArray = null;
		HashSet<LTVertex> seeds = new HashSet<LTVertex>();
		int k = 10;
		int round = 10;
		/*double ita = 0.01;
		int ell = 4;
		int N[] = {500, 600, 700, 800, 900, 1000};
		int s[] = {2, 3};
		double pc[] = new double [15];
		for (int i = 0; i < 15; i++)
		{
			pc[i] = (i+1)*0.0625;
		}
		double pm[] = new double[24];
		for (int i = 0; i < 24; i++)
		{
			pm[i] = (i+1)*0.04;
		}*/
		for (int f = 0; f < 1; f++)
		{
			vertices = new HashMap<Integer, LTVertex>();
			g = createGraph(inputFiles[f], vertices);
			System.out.println(inputFiles[f].edgesFileName+": ");
			//seeds = simpleGreedy(g, vertices, k);
			long now = System.currentTimeMillis();
			
			//output.println("f: "+f);
			/*for(int N1:N)
				for(int s1:s)
					for (double pc1:pc)
						for (double pm1:pm)*/
			{
			int N1 = 900, s1 = 2, maxG = 780;
			double pc1 = 0.5625, pm1 = 0.12;
			
			//maxG = 30000/N1;
			//seedsArray = GA(1, g, vertices, new HashSet<LTVertex>(), new HashSet<LTVertex>(), new HashSet<LTVertex>(), new HashSet<LTVertex>(), k, round, N1, s1, pc1, pm1, maxG);
			//seeds = SIMPATH(g, vertices, ita, ell, k);

			HashSet<LTVertex> activeVertices = new HashSet<LTVertex>();
			Pair<HashSet<LTVertex>> newerActiveVertices = new Pair<HashSet<LTVertex>>(new HashSet<LTVertex>(), new HashSet<LTVertex>());
			HashSet<Vertex> activatedVertices= new HashSet<Vertex>(activeVertices.size());
			
			resetState(vertices, new HashSet<LTVertex>(), new HashSet<LTVertex>(), new HashSet<LTVertex>(), new HashSet<LTVertex>());
			//resetHold(g, vertices);
			for (int r = 0; r < round; r++)
			{
				seeds.clear();
				seeds.addAll(newerActiveVertices.getFirst());
				seeds.addAll(seedsArray.get(r));
				activeVertices.clear();
				newerActiveVertices = null;
				
				newerActiveVertices = runDiffusion(g, vertices, seeds, new HashSet<LTVertex>(), 1, activeVertices, new HashSet<LTVertex>(), Model.LTModel);
				activatedVertices.addAll(activeVertices);
				activatedVertices.addAll(newerActiveVertices.getFirst());
			}
			seeds.clear();
			seeds.addAll(newerActiveVertices.getFirst());
			activeVertices.clear();
			newerActiveVertices = null;
			
			newerActiveVertices = runDiffusion(g, vertices, seeds, new HashSet<LTVertex>(), 0, activeVertices, new HashSet<LTVertex>(), Model.LTModel);
			activatedVertices.addAll(activeVertices);
			activatedVertices.addAll(newerActiveVertices.getFirst());
				System.out.println("Spread: "+activatedVertices.size());
				//output.println(""+N1+", "+s1+", "+pc1+", "+pm1+", "+maxG+", "+activatedVertices.size()+", "+(System.currentTimeMillis()-now));
				now = System.currentTimeMillis();
				//System.out.println();
			}
			
		}
		//output.close();
	}
	
	private static void testStratergyForOnePartyGreedy()
	{
		InputFileCollect[] inputFiles = new InputFileCollect[]{
				new InputFileCollect("partB_hepth_lt_edges.txt", "partB_hepth_lt_nodes.txt", ""),
				new InputFileCollect("partB_egofb_lt_edges.txt", "partB_egofb_lt_nodes.txt", ""),
		};
		
		Graph<LTVertex, LTEdge> g;
		HashMap<Integer, LTVertex> vertices;
		HashSet<LTVertex> seeds = new HashSet<LTVertex>();
		int k = 10;
		int round = 10;
		for (int f = 0; f < 2; f++)
		{
			vertices = new HashMap<Integer, LTVertex>();
			g = createGraph(inputFiles[f], vertices);
			System.out.println(inputFiles[f].edgesFileName+": ");

			HashSet<LTVertex> activeVertices = new HashSet<LTVertex>();
			Pair<HashSet<LTVertex>> newerActiveVertices = new Pair<HashSet<LTVertex>>(new HashSet<LTVertex>(), new HashSet<LTVertex>());
			HashSet<LTVertex> activatedVertices= new HashSet<LTVertex>(activeVertices.size());
			
			
			for (int r = 0; r < round; r++)
			{
				//seeds = simpleGreedy(g, vertices, k, newerActiveVertices, activatedVertices);
				resetState(vertices, newerActiveVertices.getFirst(), newerActiveVertices.getSecond(), activatedVertices, new HashSet<LTVertex>());
				//resetHold(g, vertices);
				activatedVertices.addAll(newerActiveVertices.getFirst());
				seeds.addAll(newerActiveVertices.getFirst());
				activeVertices.clear();
				newerActiveVertices = null;
				
				newerActiveVertices = runDiffusion(g, vertices, seeds, new HashSet<LTVertex>(), 1, activeVertices, new HashSet<LTVertex>(), Model.LTModel);
				activatedVertices.addAll(activeVertices);
			}
			activatedVertices.addAll(newerActiveVertices.getFirst());
			seeds.clear();
			seeds.addAll(newerActiveVertices.getFirst());
			activeVertices.clear();
			newerActiveVertices = null;
			
			newerActiveVertices = runDiffusion(g, vertices, seeds, new HashSet<LTVertex>(), 0, activeVertices, new HashSet<LTVertex>(), Model.LTModel);
			activatedVertices.addAll(activeVertices);
			activatedVertices.addAll(newerActiveVertices.getFirst());
				System.out.println("Spread: "+activatedVertices.size());
				//output.println(""+N1+", "+s1+", "+pc1+", "+pm1+", "+maxG+", "+activatedVertices.size()+", "+(System.currentTimeMillis()-now));
				//System.out.println();
		}
	}
	
	private static ArrayList<HashSet<LTVertex>> GA(int playerID, Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, HashMap<LTVertex, Pair<Double>> ltHoldMap,
			HashSet<LTVertex> player1NewlyActiveVertices, HashSet<LTVertex> player2NewlyActiveVertices, 
			HashSet<LTVertex> player1ActiveVertices, HashSet<LTVertex> player2ActiveVertices,
			int k, int round, int N, int s, double pc, double pm, int timeLimit)
	{
		int saveChromosomeNum = 5;
		ArrayList<HashSet<LTVertex>> S = new ArrayList<HashSet<LTVertex>>(saveChromosomeNum+1);
		for (int r = 0; r < saveChromosomeNum+1; r++)
			S.add(new HashSet<LTVertex>(k));
		SNA_FitnessFunction f = new SNA_FitnessFunction(playerID, g, vertices, ltHoldMap, 1,
				player1NewlyActiveVertices, player2NewlyActiveVertices, player1ActiveVertices, player2ActiveVertices);
		
		GA ga = new GA(GA.SelectionModel.TOURNAMENT_SELECT, f, vertices.size(), N, s, pc, pm, timeLimit, -1, k, 1, vertices.values().toArray(new LTVertex[0]), saveChromosomeNum, round == 10);
		ga.doIt(false);
		ga.showStatistics();
		
		for (int i = 0; i < ga.nInitial; i++)
			ga.selectionIndex[i] = i;
		quickSort(ga.selectionIndex, ga.population, 0, ga.nInitial-1);
		 
		Chromosome pChromo = null;
		int count = 0;
		int countSaved = 0;
		for (countSaved = 0; countSaved < saveChromosomeNum+1; countSaved++)
		{
			while (count < 1000 && pChromo != null && pChromo.gene.get(0).containsAll(ga.population[ga.selectionIndex[ga.nInitial-1-count]].gene.get(0)))
				count++;
			if (count == 1000)
				break;
			pChromo = ga.population[ga.selectionIndex[ga.nInitial-1-count]];
			count++;
			for (Integer ind:pChromo.gene.get(0))
			{
				S.get(countSaved).add(pChromo.arrVertices[ind]);
			}
		}
		PrintStream output = null;
		try {
			output = new PrintStream(new File("GA_preSelectedGene.txt"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (int r = 1; r < saveChromosomeNum && r < countSaved; r++)
		{
			for (LTVertex v:S.get(r))
				output.print(v.vertexIndex+" ");
			output.println();
		}
		output.close();
		return S;
	}
	
	/*private static HashSet<LTVertex> simpleGreedy(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, int k, Pair<HashSet<LTVertex>> origNewerActiveVertices,
			HashSet<LTVertex> origActivatedVertices)
	{
		HashSet<LTVertex> S = new HashSet<LTVertex>(k);
		List<VertexWithData> CELF_queue = new SortedVertexList();
		HashSet<LTVertex> activeVertices = new HashSet<LTVertex>(g.getVertexCount());
		HashSet<LTVertex> activatedVertices = null;
		Pair<HashSet<LTVertex>> newerActiveVertices = null;
		LTVertex best;
		VertexWithData bestWithData;
		
		int spread = 0;
		for (LTVertex v:vertices.values())
		{
			S.add(v);
			activeVertices.clear();
			resetState(vertices, origNewerActiveVertices.getFirst(), origNewerActiveVertices.getSecond(), origActivatedVertices, new HashSet<LTVertex>());
			resetHold(g, vertices);
			newerActiveVertices = runDiffusion(g, vertices, S, new HashSet<LTVertex>(), 0, activeVertices, new HashSet<LTVertex>(), Model.LTModel);
			activatedVertices = new HashSet<LTVertex>(activeVertices.size());
			activatedVertices.addAll(activeVertices);
			activatedVertices.addAll(newerActiveVertices.getFirst());
			CELF_queue.add(new VertexWithData(v, activatedVertices.size()-spread));
			S.remove(v);
		}
		while (S.size()<k)
		{
			best = CELF_queue.get(CELF_queue.size()-1).vertex;
			CELF_queue.remove(CELF_queue.size()-1);
			S.add(best);
			resetState(vertices, origNewerActiveVertices.getFirst(), origNewerActiveVertices.getSecond(), origActivatedVertices, new HashSet<LTVertex>());
			resetHold(g, vertices);
			activeVertices.clear();
			newerActiveVertices = runDiffusion(g, vertices, S, new HashSet<LTVertex>(), 0, activeVertices, new HashSet<LTVertex>(), Model.LTModel);
			S.remove(best);
			activatedVertices.clear();
			activatedVertices.addAll(newerActiveVertices.getFirst());
			activatedVertices.addAll(activeVertices);
			bestWithData = new VertexWithData(best, activatedVertices.size()-spread);
			CELF_queue.add(bestWithData);
			if (CELF_queue.get(CELF_queue.size()-1).data == bestWithData.data)
			{
				CELF_queue.remove(bestWithData);
				S.add(best);
			}
		}
		return S;
	}
	
	private static HashSet<LTVertex> SIMPATH(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, double ita, int ell, int k)
	{
		HashSet<LTVertex> C = findVC(g, vertices);
		HashSet<LTVertex> U = new HashSet<LTVertex>();
		HashSet<LTVertex> S = new HashSet<LTVertex>(k);
		HashMap<LTVertex, Double> spreadMinusV = new HashMap<LTVertex, Double>();
		HashMap<LTVertex, HashMap<LTVertex, Double>> spreadMinusV_of_u = new HashMap<LTVertex, HashMap<LTVertex, Double>>(vertices.size());
		HashMap<LTVertex, Double> marginalGain = new HashMap<LTVertex, Double>(vertices.size());
		LTVertex arrVertices[] = vertices.values().toArray(new LTVertex[0]);
		double arrMarginalGain[] = new double[vertices.size()];
		int sortedIndex[] = new int[vertices.size()];
		for (int i = 0; i < vertices.size(); i++)
		{
			sortedIndex[i] = i;
		}
		
		HashSet<LTVertex> V_minus_C = new HashSet<LTVertex>();
		V_minus_C.addAll(vertices.values());
		V_minus_C.removeAll(C);
		for (LTVertex u:C)
		{
			U.clear();
			U.addAll(V_minus_C);
			U.retainAll(g.getPredecessors(u));
			S.clear();
			S.add(u);
			spreadMinusV.clear();
			marginalGain.put(u, SIMPATH_spread(g, vertices.values(), S, ita, U, spreadMinusV));
			spreadMinusV_of_u.put(u, spreadMinusV);
		}
		
		
		double spreadV;
		for (LTVertex v:V_minus_C)
		{
			spreadV = 1;
			for (LTVertex u:g.getSuccessors(v))
			{
				if (spreadMinusV_of_u.get(u).containsKey(v))
					spreadV += g.findEdge(v, u).ltInfluence*spreadMinusV_of_u.get(u).get(v);
			}
			marginalGain.put(v, spreadV);
		}
		
		S.clear();
		double spread = 0;
		HashSet<LTVertex> V_minus_S = new HashSet<LTVertex>();
		HashSet<LTVertex> examined = new HashSet<LTVertex>();
		double spread_V_minus_S_of_x;
		while (S.size() < k)
		{
			for (int i = 0; i < arrVertices.length; i++)
				arrMarginalGain[i] = marginalGain.get(arrVertices[i]);
			quickSort(sortedIndex, arrMarginalGain, 0, arrMarginalGain.length - 1);
			U.clear();
			for (int i = 0; i < ell; i++)
				U.add(arrVertices[sortedIndex[arrVertices.length-1-i]]);
			spreadMinusV.clear();
			SIMPATH_spread(g, vertices.values(), S, ita, U, spreadMinusV);
			for (LTVertex x:U)
			{
				if (examined.contains(x))
				{
					S.add(x);
					spread += marginalGain.get(x);
					examined.clear();
				}
				V_minus_S.clear();
				V_minus_S.addAll(vertices.values());
				V_minus_S.removeAll(S);
				spread_V_minus_S_of_x = backtrack(g, x, ita, V_minus_S, null, null);
				if (!spreadMinusV.containsKey(x))
					marginalGain.replace(x, spread_V_minus_S_of_x - spread);
				else
				{
					marginalGain.replace(x, spreadMinusV.get(x)+spread_V_minus_S_of_x - spread);
				}
				examined.add(x);
				
			}
		}
			
		
		
		
		return S;
	}
	
	private static HashSet<LTVertex> findVC(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices)
	{
		HashSet<LTVertex> C = new HashSet<LTVertex>();
		List<VertexWithData> heuristicWithDegree = new SortedVertexList();
		HashMap<Vertex, VertexWithData> map = new HashMap<Vertex, VertexWithData>();
		VertexWithData tempVwithD;
		HashSet<LTEdge> edges = new HashSet<LTEdge>();
		LTVertex vertexToC;
		Collection<LTEdge> tempIncidentEdges;
		edges.addAll(g.getEdges());
		
		for (LTVertex v:vertices.values())
		{
			tempVwithD = new VertexWithData(v, g.degree(v));
			heuristicWithDegree.add(tempVwithD);
			map.put(v, tempVwithD);
		}
		while (!edges.isEmpty())
		{
			vertexToC = heuristicWithDegree.get(heuristicWithDegree.size()-1).vertex;
			tempIncidentEdges = g.getIncidentEdges(vertexToC);
			if (edges.removeAll(tempIncidentEdges))			//remove covered edges. if succeed, vertexToC is chose.
			{
				C.add(vertexToC);
				for (LTEdge e:tempIncidentEdges)
				{
					heuristicWithDegree.remove(map.get(g.getOpposite(vertexToC, e)));
					map.get(g.getOpposite(vertexToC, e)).data--;
					heuristicWithDegree.add(map.get(g.getOpposite(vertexToC, e)));
				}
				heuristicWithDegree.remove(map.get(vertexToC));
			}
			else
			{
				heuristicWithDegree.remove(map.get(vertexToC));
			}
		}
		return C;
	}
	
	private static double SIMPATH_spread(Graph<LTVertex, LTEdge> g, Collection<LTVertex> vertices, HashSet<LTVertex> seeds, double ita,
			HashSet<LTVertex> U, HashMap<LTVertex, Double> spreadMinusV)
	{
		int spread = 0;
		HashSet<LTVertex> W = new HashSet<LTVertex>(vertices.size()-seeds.size()+1);
		W.addAll(vertices);
		W.removeAll(seeds);
		for (LTVertex u:seeds)
		{
			//W = V - S + u
			W.add(u);
			spread += backtrack(g, u, ita, W, U, spreadMinusV);
			W.remove(u);
		}
		return spread;
	}
	
	private static double backtrack(Graph<LTVertex, LTEdge> g, LTVertex u, double ita,
			HashSet<LTVertex> W, HashSet<LTVertex> U, HashMap<LTVertex, Double> spreadMinusV)
	{
		Stack<LTVertex> Q = new Stack<LTVertex>();
		DoubleRef spread = new DoubleRef(1);
		DoubleRef pp = new DoubleRef(1);
		HashMap<LTVertex, HashSet<LTVertex>> D = new HashMap<LTVertex, HashSet<LTVertex>>();
		LTVertex v;
		Q.push(u);
		
		while(!Q.empty())
		{
			forward(g, Q, D, spread, pp, ita, W, U, spreadMinusV);
			u = Q.pop();
			D.remove(u);
			if (Q.empty())
				break;
			v = Q.peek();
			pp.value /= g.findEdge(v, u).ltInfluence;
		}
		
		return spread.value;
	}

	private static void forward(Graph<LTVertex, LTEdge> g, Stack<LTVertex> Q, HashMap<LTVertex, HashSet<LTVertex>> D, 
		DoubleRef spread, DoubleRef pp, double ita, HashSet<LTVertex> W, HashSet<LTVertex> U, HashMap<LTVertex, Double> spreadMinusV)
	{
		LTVertex x = Q.peek();
		LTVertex y;
		HashSet<LTVertex> outNeighbor = new HashSet<LTVertex>();
		
		outNeighbor.addAll(g.getSuccessors(x));
		outNeighbor.removeAll(Q);
		if (D.containsKey(x))
			outNeighbor.removeAll(D.get(x));
		outNeighbor.retainAll(W);
		while(!outNeighbor.isEmpty())
		{
			y = outNeighbor.iterator().next();
			if (pp.value*g.findEdge(x, y).ltInfluence < ita)
			{
				insert(D, x, y);
				outNeighbor.remove(y);
			}
			else
			{
				Q.add(y);
				pp.value *= g.findEdge(x, y).ltInfluence;
				spread.value += pp.value;
				insert(D, x, y);
				x = Q.peek();
				if (U != null)
				{
					for (LTVertex v:U)
					{
						if (!Q.contains(v))
						{
							if (!spreadMinusV.containsKey(v))
								spreadMinusV.put(v, pp.value);
							else
								spreadMinusV.replace(v, spreadMinusV.get(v) + pp.value);
						}
					}
				}
				outNeighbor.clear();
				outNeighbor.addAll(g.getSuccessors(x));
				outNeighbor.removeAll(Q);
				if (D.containsKey(x))
					outNeighbor.removeAll(D.get(x));
				outNeighbor.retainAll(W);
			}
		}
	}
	
	private static void insert(HashMap<LTVertex, HashSet<LTVertex>> D, LTVertex x, LTVertex y)
	{
		if (!D.containsKey(x))
			D.put(x, new HashSet<LTVertex>());
		D.get(x).add(y);
	}
	*/
	private static Graph<LTVertex, LTEdge> createGraph(InputFileCollect inputFile, HashMap<Integer, LTVertex> vertices)
	{
		BufferedReader input = null;
		String line = null;
		String nodes[] = new String[3];
		LTVertex vertex1, vertex2;
		Graph<LTVertex, LTEdge> g = null;
		//Create Graph g
		try 
		{
			input = new BufferedReader(new FileReader(inputFile.edgesFileName));
			line = input.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (line.equals("# Directed"))
			g = new DirectedSparseMultigraph<LTVertex, LTEdge>();
		else if (line.equals("# Undirected"))
			g = new UndirectedSparseMultigraph<LTVertex, LTEdge>();
		else
			System.out.println("Wrong input file!!");
		try 
		{
			input.readLine();
			line = input.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		while (line != null)
		{
			nodes = line.split(" ");
			vertex1 = getVertex(Integer.parseInt(nodes[0]), vertices, inputFile.model);
			addVertex(vertex1, g, vertices);
			vertex2 = getVertex(Integer.parseInt(nodes[1]), vertices, inputFile.model);
			addVertex(vertex2, g, vertices);
				g.addEdge(new LTEdge(Double.parseDouble(nodes[2])), vertex1, vertex2);
			
			try 
			{
				line = input.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		try {
			input.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if (inputFile.model == Model.LTModel)
		{
			String term[];
			try 
			{
				input = new BufferedReader(new FileReader(inputFile.nodesFileName));
				input.readLine();
				line = input.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			while (line != null)
			{
				term = line.split(" ");
				((LTVertex)vertices.get(Integer.parseInt(term[0]))).ltThreshold = Double.parseDouble(term[1]);
				try {
					line = input.readLine();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			try {
				input.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return g;
	}
	
	private static void resetState(HashMap<Integer, LTVertex> vertices, 
			HashSet<LTVertex> player1NewlyActiveVertices, HashSet<LTVertex> player2NewlyActiveVertices, 
			HashSet<LTVertex> player1ActiveVertices, HashSet<LTVertex> player2ActiveVertices)
	{
		for (Vertex v:vertices.values())
		{
			if (player1NewlyActiveVertices.contains(v))
				v.state = Vertex.State.PLAYER1_NEWLY_ACTIVE;
			else if (player2NewlyActiveVertices.contains(v))
				v.state = Vertex.State.PLAYER2_NEWLY_ACTIVE;
			else if (player1ActiveVertices.contains(v))
				v.state = Vertex.State.PLAYER1_ACTIVE;
			else if (player2ActiveVertices.contains(v))
				v.state = Vertex.State.PLAYER2_ACTIVE;
			else
				v.state = Vertex.State.INACTIVE;
		}
	}
	
	private static void resetHold(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, HashMap<LTVertex, Pair<Double>> ltHoldMap)
	{
		for (LTVertex v:vertices.values())
		{
			if (ltHoldMap != null)
			{
				v.ltHold1 = ltHoldMap.get(v).getFirst();
				v.ltHold2 = ltHoldMap.get(v).getSecond();
			}
			else
			{
				v.ltHold1 = 0;
				v.ltHold2 = 0;
			}
		}
	}
	
	private static HashMap<LTVertex, Pair<Double>> createLTHoldMap(HashMap<Integer, LTVertex> vertices)
	{
		HashMap<LTVertex, Pair<Double>> holdMap = new HashMap<LTVertex, Pair<Double>>(vertices.size());
		BufferedReader input = null;
		String line = null;
		String token[];
		try 
		{
			input = new BufferedReader(new FileReader("GA_holdStatus.txt"));
			line = input.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		while(line != null)
		{
			token = line.split(" ");
			holdMap.put(vertices.get(Integer.parseInt(token[0])), new Pair<Double>(Double.parseDouble(token[1]), Double.parseDouble(token[2])));
			try {
				line = input.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		try {
			input.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return holdMap;
	}
	
	private static void updateLTHold(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, 
			HashSet<LTVertex> player1LastChoseVertices, HashSet<LTVertex> player2LastChoseVertices)
	{
		PrintStream output = null;
		try {
			output = new PrintStream(new File("GA_holdStatus.txt"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		HashSet<LTVertex> activeVertices1 = new HashSet<LTVertex>();
		HashSet<LTVertex> activeVertices2 = new HashSet<LTVertex>();

		runDiffusion(g, vertices, player1LastChoseVertices, player2LastChoseVertices, 0, activeVertices1, activeVertices2, Model.LTModel);
		
		for (LTVertex v:vertices.values())
			output.println(v.vertexIndex+" "+v.ltHold1+" "+v.ltHold2);
		output.close();
	}
	
	private static LTVertex getVertex(int vertexIndex, HashMap<Integer, LTVertex> vertices, Model model)
	{
		if (vertices.containsKey(vertexIndex))
			return vertices.get(vertexIndex);
		else
		{
				return new LTVertex(vertexIndex, 0, Vertex.State.INACTIVE);
		}
	}
	
	private static void addVertex(LTVertex v, Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices)
	{
		if (!vertices.containsKey(v.vertexIndex))
		{
			g.addVertex(v);
			vertices.put(v.vertexIndex, v);
		}
	}
	
	private static int partition(int index[], double arr[], int left, int right)
	{
	      int i = left, j = right;
	      int tmp;
	      double pivot = arr[index[(left + right) / 2]];
	     
	      while (i <= j) {
	            while (arr[index[i]] < pivot)
	                  i++;
	            while (arr[index[j]] > pivot)
	                  j--;
	            if (i <= j) {
	                  tmp = index[i];
	                  index[i] = index[j];
	                  index[j] = tmp;
	                  i++;
	                  j--;
	            }
	      };
	     
	      return i;
	}
	 
	private static void quickSort(int index[], double arr[], int left, int right) {
	      int indexI = partition(index, arr, left, right);
	      if (left < indexI - 1)
	            quickSort(index, arr, left, indexI - 1);
	      if (indexI < right)
	            quickSort(index, arr, indexI, right);
	}
	
	private static int partition(int index[], Chromosome arr[], int left, int right)
	{
	      int i = left, j = right;
	      int tmp;
	      Chromosome pivot = arr[index[(left + right) / 2]];
	     
	      while (i <= j) {
	            while (arr[index[i]].getFitness() < pivot.getFitness())
	                  i++;
	            while (arr[index[j]].getFitness() > pivot.getFitness())
	                  j--;
	            if (i <= j) {
	                  tmp = index[i];
	                  index[i] = index[j];
	                  index[j] = tmp;
	                  i++;
	                  j--;
	            }
	      };
	     
	      return i;
	}
	 
	private static void quickSort(int index[], Chromosome arr[], int left, int right) {
	      int indexI = partition(index, arr, left, right);
	      if (left < indexI - 1)
	            quickSort(index, arr, left, indexI - 1);
	      if (indexI < right)
	            quickSort(index, arr, indexI, right);
	}
	
	private static Pair<HashSet<LTVertex>> runDiffusion(Graph<LTVertex, LTEdge> g, HashMap<Integer, LTVertex> vertices, 
			HashSet<LTVertex> newlyActiveVertices1, HashSet<LTVertex> newlyActiveVertices2, int rounds, 
			HashSet<LTVertex> activeVertices1, HashSet<LTVertex> activeVertices2, Model model)
	{
		// TODO Auto-generated method stub
		Collection<LTEdge> edges;
		LTVertex successor;
		HashSet<LTVertex> newerActiveVertices1 = null;
		HashSet<LTVertex> newerActiveVertices2 = null;
		for (int round = 0; round < rounds || (rounds==0); round++)
		{
			newerActiveVertices1 = new HashSet<LTVertex>(g.getVertexCount());
			newerActiveVertices2 = new HashSet<LTVertex>(g.getVertexCount());
			for (LTVertex v:newlyActiveVertices1)
			{
				if (v.state == Vertex.State.INACTIVE||v.state == Vertex.State.PLAYER1_NEWLY_ACTIVE)
				{
					v.state = Vertex.State.PLAYER1_ACTIVE;
					activeVertices1.add(v);
					edges = g.getOutEdges(v);
					for (LTEdge edge:edges)
					{
						successor = g.getDest(edge);
						if ((successor.state == Vertex.State.PLAYER1_NEWLY_ACTIVE || successor.state == Vertex.State.INACTIVE) && !newlyActiveVertices2.contains(successor))
						{
							((LTVertex)successor).ltHold1 += ((LTEdge)edge).ltInfluence;
							if (((LTVertex)successor).ltHold1 >= ((LTVertex)successor).ltThreshold)
							{
								successor.state = Vertex.State.PLAYER1_NEWLY_ACTIVE;
								newerActiveVertices1.add(successor);
							}
						}
					}
				}
			}
			for (LTVertex v:newlyActiveVertices2)
			{
				if (v.state == Vertex.State.INACTIVE||v.state == Vertex.State.PLAYER2_NEWLY_ACTIVE)
				{
					v.state = Vertex.State.PLAYER2_ACTIVE;
					activeVertices2.add(v);
					edges = g.getOutEdges(v);
					for (LTEdge edge:edges)
					{
						successor = g.getDest(edge);
						if (successor.state == Vertex.State.PLAYER2_NEWLY_ACTIVE || successor.state == Vertex.State.INACTIVE)
						{
							((LTVertex)successor).ltHold2 += ((LTEdge)edge).ltInfluence;
							if (((LTVertex)successor).ltHold2 >= ((LTVertex)successor).ltThreshold)
							{
								successor.state = Vertex.State.PLAYER2_NEWLY_ACTIVE;
								newerActiveVertices2.add(successor);
							}
						}
						else if (successor.state == Vertex.State.PLAYER1_NEWLY_ACTIVE)
						{
							((LTVertex)successor).ltHold2 += ((LTEdge)edge).ltInfluence;
							if (((LTVertex)successor).ltHold2 >= ((LTVertex)successor).ltHold1 
									&& ((LTVertex)successor).ltHold2 >= ((LTVertex)successor).ltThreshold)
							{
								successor.state = Vertex.State.PLAYER2_NEWLY_ACTIVE;
								newerActiveVertices2.add(successor);
								newerActiveVertices1.remove(successor);
							}
						}
					}
				}
			}
			if (newerActiveVertices1.size() == 0 && newerActiveVertices2.size() == 0)
			{
				//System.out.println("End of diffusion when round"+round);
				break;
			}
			newlyActiveVertices1 = newerActiveVertices1;
			newlyActiveVertices2 = newerActiveVertices2;
		}
		
		return new Pair<HashSet<LTVertex>>(newerActiveVertices1, newerActiveVertices2);
	}
}