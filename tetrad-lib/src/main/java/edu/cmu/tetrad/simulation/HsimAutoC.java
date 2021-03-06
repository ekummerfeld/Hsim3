package edu.cmu.tetrad.simulation;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.io.TabularContinuousDataReader;
import edu.cmu.tetrad.io.VerticalTabularDiscreteDataReader;
import edu.cmu.tetrad.search.*;

import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by Erich on 3/28/2016.
 */
public class HsimAutoC {
    private boolean verbose=false;
    private DataSet data;
    private boolean write = false;
    private String filenameOut = "defaultOut";
    private char delimiter = ',';

    //*********Constructors*************//
    //contructor using a previously existing DataSet object
    public HsimAutoC(DataSet indata) {
        //first check if indata is already the right type
        data = indata;
        //may need to make this part more complicated if CovarianceMatrix method is finicky
    }

    //constructor that loads data from a file named readfilename, with delimiter delim
    public HsimAutoC(String readfilename, char delim){
        String workingDirectory = System.getProperty("user.dir");
        System.out.println(workingDirectory);
        Set<String> eVars = new HashSet<String>();
        eVars.add("MULT");
        Path dataFile = Paths.get(readfilename);

        edu.cmu.tetrad.io.DataReader dataReader = new TabularContinuousDataReader(dataFile, delim);
        try {
            data = dataReader.readInData(eVars);
        }
        catch(Exception IOException){
            IOException.printStackTrace();
        }
    }

    //***********Public methods*************//

    public double[] run(int resimSize) {

        double[] output;
        output = new double[5];
        //========first make the Dag for Hsim==========
        ICovarianceMatrix cov = new CovarianceMatrixOnTheFly(data);
        SemBicScore score = new SemBicScore(cov);

        double penaltyDiscount = 2.0;
        Fges fges = new Fges(score);
        fges.setVerbose(false);
        fges.setNumPatternsToStore(0);
        fges.setPenaltyDiscount(penaltyDiscount);

        Graph estGraph = fges.search();
        //if (verbose) System.out.println(estGraph);
        estGraph = dagFromPattern(estGraph);
        Dag estDAG = new Dag(estGraph);
        //Dag estDAG = new Dag(estGraph);

        //===========Identify the nodes to be resimulated===========

        //for this class, I'm going to choose variables for resimulation randomly, rather than building cliques

        //select a random node
        List<Node> remainingNodes = estGraph.getNodes();
        int randIndex = new Random().nextInt(remainingNodes.size());
        Node randomnode = remainingNodes.get(randIndex);
        if (verbose) System.out.println("the first node is " + randomnode);
        List<Node> queue = new ArrayList<>();
        queue.add(randomnode);
        //while queue has size less than the resim size, grow it
        //if (verbose) System.out.println(queue);
        while (queue.size() < resimSize) {
            //choose another node randomly
            remainingNodes.remove(randIndex);
            randIndex = new Random().nextInt(remainingNodes.size());
            randomnode = remainingNodes.get(randIndex);
            //add that node to the resim set
            queue.add(randomnode);
        }

        Set<Node> simnodes = new HashSet<Node>(queue);
        if (verbose) System.out.println("the resimmed nodes are " + simnodes);

        //===========Apply the hybrid resimulation===============
        HsimContinuous hsimC = new HsimContinuous(estDAG,simnodes,data); //regularDataSet
        DataSet newDataSet = hsimC.hybridsimulate();

        //write output to a new file

        if (write) {
            try {
                FileWriter fileWriter = new FileWriter(filenameOut);
                DataWriter.writeRectangularData(newDataSet, fileWriter, delimiter);
                fileWriter.close();
            }
            catch(Exception IOException){
                IOException.printStackTrace();
            }
        }

        //=======Run FGS on the output data, and compare it to the original learned graph
        //Path dataFileOut = Paths.get(filenameOut);
        //edu.cmu.tetrad.io.DataReader dataReaderOut = new VerticalTabularDiscreteDataReader(dataFileOut, delimiter);

        ICovarianceMatrix newcov = new CovarianceMatrixOnTheFly(data);
        SemBicScore newscore = new SemBicScore(newcov);
        Fges fgesOut = new Fges(newscore);
        fgesOut.setVerbose(false);
        fgesOut.setNumPatternsToStore(0);
        fgesOut.setPenaltyDiscount(2.0);

        Graph estGraphOut = fgesOut.search();
        //if (verbose) System.out.println(" bugchecking: fgs estGraphOut: " + estGraphOut);

        //doing the replaceNodes trick to fix some bugs
        estGraphOut = GraphUtils.replaceNodes(estGraphOut,estDAG.getNodes());
        //restrict the comparison to the simnodes and edges to their parents
        Set<Node> allParents = HsimUtils.getAllParents(estGraphOut, simnodes);
        Set<Node> addParents = HsimUtils.getAllParents(estDAG, simnodes);
        allParents.addAll(addParents);
        Graph estEvalGraphOut = HsimUtils.evalEdges(estGraphOut,simnodes,allParents);
        Graph estEvalGraph = HsimUtils.evalEdges(estDAG,simnodes,allParents);

        //SearchGraphUtils.graphComparison(estGraph, estGraphOut, System.out);

        estEvalGraphOut = GraphUtils.replaceNodes(estEvalGraphOut, estEvalGraph.getNodes());
        //if (verbose) System.out.println(estEvalGraph);
        //if (verbose) System.out.println(estEvalGraphOut);

        //SearchGraphUtils.graphComparison(estEvalGraphOut, estEvalGraph, System.out);
        output = HsimUtils.errorEval(estEvalGraphOut,estEvalGraph);
        if (verbose) System.out.println(output[0]+" "+output[1]+" "+output[2]+" "+output[3]+" "+output[4]);
        return output;
    }
    private static Graph dagFromPattern(Graph pattern) {
        Graph dag = new EdgeListGraph(pattern);;
        MeekRules rules = new MeekRules();

        Random rand = new Random();
        WHILE:
        while (true) {
            List<Edge> edges = new ArrayList<>(dag.getEdges());
            Collections.shuffle(edges);
            for (Edge edge : edges) {
                if (Edges.isUndirectedEdge(edge)) {
                    Node x, y;
                    if (rand.nextBoolean()){
                        x = edge.getNode1();
                        y = edge.getNode2();
                    } else {
                        y = edge.getNode1();
                        x = edge.getNode2();
                    }

                    List<Node> okx = dag.getAdjacentNodes(x);
                    okx.removeAll(dag.getChildren(x));
                    okx.remove(y);

                    List<Node> oky = dag.getAdjacentNodes(y);
                    oky.removeAll(dag.getChildren(y));
                    oky.remove(x);

                    if (!okx.isEmpty()) {
                        Node other = okx.get(0);
                        dag.removeEdge(other, x);
                        dag.removeEdge(y, x);
                        dag.addDirectedEdge(other, x);
                        dag.addDirectedEdge(y, x);
                    } else if (!oky.isEmpty()) {
                        Node other = oky.get(0);
                        dag.removeEdge(other, y);
                        dag.removeEdge(x, y);
                        dag.addDirectedEdge(other, y);
                        dag.addDirectedEdge(x, y);
                    } else {
                        dag.removeEdge(x, y);
                        dag.addDirectedEdge(x, y);
                    }

                    rules.orientImplied(dag);
                    continue WHILE;
                }
            }
            if (!dag.existsDirectedCycle()){
                break;
            } else {
                dag = new EdgeListGraph(pattern);
            }
        }

        return dag;
    }
    //******* Methods for setting values to private variables****************//
    public void setVerbose(boolean verbosity){
        verbose=verbosity;
    }
    public void setWrite(boolean setwrite) {write=setwrite;}
    public void setFilenameOut(String filename) { filenameOut=filename;}
    public void setDelimiter(char delim){delimiter=delim;}
}
