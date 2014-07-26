package tips.dc;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import com.beust.jcommander.internal.Lists;

import briefj.collections.Counter;
import briefj.collections.Tree;
import briefj.collections.UnorderedPair;
import conifer.TreeNode;
import conifer.UnrootedTree;



public class TipsTreeUtils
{
  public static Tree<Pair<TreeNode,Double>> topologicalCentroidRooting(UnrootedTree tree)
  {
    // find a centroid
    Counter<TreeNode> internalTipsDistances = internalTipsDistances(tree);
    TreeNode multifurcatingCentroid = internalTipsDistances.argMin();
    
    // find a min neighbor
    double min = Double.POSITIVE_INFINITY;
    TreeNode argmin = null;
    for (TreeNode neighbor : Graphs.neighborListOf(tree.getTopology(), multifurcatingCentroid))
      if (argmin == null || internalTipsDistances.getCount(neighbor) < min)
      {
        argmin = neighbor;
        min = internalTipsDistances.getCount(neighbor);
      }
    
    double halfBL = tree.getBranchLength(multifurcatingCentroid, argmin) / 2.0;
    Tree<Pair<TreeNode,Double>> 
      left = topologicalCentroidRooting(tree, multifurcatingCentroid, argmin, halfBL),
      right= topologicalCentroidRooting(tree, argmin, multifurcatingCentroid, halfBL);
    
    @SuppressWarnings("unchecked")
    List<Tree<Pair<TreeNode,Double>>> children = Arrays.asList(left, right);
    Pair<TreeNode, Double> label = Pair.of(TreeNode.nextUnlabelled(),null);
    return new Tree<Pair<TreeNode,Double>>(label, children);
  }
 
  
  private static Tree<Pair<TreeNode, Double>> topologicalCentroidRooting(
      UnrootedTree tree, TreeNode node, TreeNode parent, double branchLen)
  {
    Pair<TreeNode, Double> pair = Pair.of(node, branchLen);
    List<Tree<Pair<TreeNode,Double>>> children = Lists.newArrayList();
    
    for (TreeNode neighbor : Graphs.neighborListOf(tree.getTopology(), node))
      if (!neighbor.equals(parent))
        children.add(topologicalCentroidRooting(tree, neighbor, node, tree.getBranchLength(neighbor, node)));
    
    return new Tree<Pair<TreeNode,Double>>(pair, children);
  }


  /**
   * Return for each node sum of the number of edge to each leaves.
   * @param tree
   * @return
   */
  public static Counter<TreeNode> internalTipsDistances(UnrootedTree tree)
  {
    Counter<TreeNode> result = new Counter<TreeNode>();
    
    for (TreeNode leaf : tree.leaves())
      addInternalTipsDistances(tree.getTopology(), leaf, null, result, 0);
    
    return result;
  }

  private static void addInternalTipsDistances(
      UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> topology, 
      TreeNode node, TreeNode parent,
      Counter<TreeNode> result, int distance)
  {
    if (parent != null)
      result.incrementCount(node, distance);
    
    for (TreeNode neighbor : Graphs.neighborListOf(topology, node))
      if (!neighbor.equals(parent))
        addInternalTipsDistances(topology, neighbor, node, result, distance + 1);
  }
  

}
