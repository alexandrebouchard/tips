package tips.dc;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


import com.google.common.collect.Lists;

import conifer.TreeNode;



public class TipsTreeObservation<S>
{
  private final List<Map<TreeNode,S>> data = Lists.newArrayList();
  
  public TipsTreeObservation(int nSites)
  {
    for (int site = 0; site < nSites; site++)
      data.add(new LinkedHashMap<TreeNode, S>());
  }

  public int nSites()
  {
    return data.size();
  }

  public S get(int siteIndex, TreeNode node)
  {
    return data.get(siteIndex).get(node);
  }
  
  public void set(int siteIndex, TreeNode node, S value)
  {
    data.get(siteIndex).put(node, value);
  }
  
  @Override
  public String toString()
  {
    return data.toString();
  }
  
}