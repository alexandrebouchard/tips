package tips.dc;


import org.junit.Assert;
import org.junit.Test;

import briefj.collections.Counter;

import conifer.TreeNode;
import conifer.UnrootedTree;



public class TestTipsTreeUtils
{
  @Test
  public void testCentroidRooting()
  {
    String treeStr =   "((A:1,B:2):3,(C:4,D:5):6);";
    UnrootedTree tree = UnrootedTree.fromNewickString(treeStr);
    tree.simplify();
    Assert.assertEquals("((unlabelled_7,null) ((unlabelled_1,4.5) (A,1.0) (B,2.0)) ((unlabelled_4,4.5) (C,4.0) (D,5.0)))".replaceAll("unlabelled[_][0-9]+", ""), TipsTreeUtils.topologicalCentroidRooting(tree).toString().replaceAll("unlabelled[_][0-9]+", ""));
  }
  
  @Test
  public void testDistance()
  {
    String treeStr =   "((A:1,(X:1,E:1):1):1,(C:1,B:2):1);";
    UnrootedTree tree = UnrootedTree.fromNewickString(treeStr);
    tree.simplify();
    System.out.println(tree.getTopology());
    System.out.println(tree);
    Counter<TreeNode> internalTipsDistances = TipsTreeUtils.internalTipsDistances(tree);
    Assert.assertEquals(internalTipsDistances.max(), 13.0, 0.0);
  }
  
  public static void main(String [] args)
  {
    new TestTipsTreeUtils().testDistance();
  }
}
