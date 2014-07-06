package tips.pip;
import java.util.ArrayList;
import java.util.List;



public class PIPString
{
  public final List<Integer> characters;
  private int [] cachedPluses = null;
  private int cachedNZeroes = Integer.MIN_VALUE;
  
  @Override
  public String toString()
  {
    return characters.toString();
  }
  
  public String toStars() 
  {
    StringBuilder result = new StringBuilder();
    for (int i = 0; i < characters.size(); i++)
      result.append("*");
    return result.toString();
  }
  
  public int [] plusses() 
  {
    ensurePlusZeroCache();
    return cachedPluses;
  }
  public int zeroes()
  {
    ensurePlusZeroCache();
    return cachedNZeroes;
  }
  
  private void ensurePlusZeroCache()
  {
    if (cachedPluses != null) return;
    cachedNZeroes = 0;
    for (int i : characters)
      if (i == 0)
        cachedNZeroes++;
    cachedPluses = new int[cachedNZeroes+1];
    int curNPlus = 0;
    int interZeroIdx = 0;
    for (int curChar : characters)
    {
      if (curChar == 0)
      {
        cachedPluses[interZeroIdx] = curNPlus;
        curNPlus = 0;
        interZeroIdx++;
      }
      else if (curChar == +1)
        curNPlus++;
      else
        throw new RuntimeException();
    }
    cachedPluses[interZeroIdx] = curNPlus;
  }
  
  
  public PIPString(List<Integer> characters)
  {
    this.characters = characters;
  }
  public PIPString(String s)
  {
    String [] split = s.split("\\s+");
    this.characters = new ArrayList<Integer>();
    for (String str : split)
      this.characters.add(Integer.parseInt(str));
  }
  
  public PIPString(List<Integer> l1, Integer pt, List<Integer> l3)
  {
    this.characters = new ArrayList<Integer>();
    this.characters.addAll(l1);
    if (pt != null) this.characters.add(pt);
    this.characters.addAll(l3);
  } 
  public PIPString(List<Integer> l1, List<Integer> l2)
  {
    this(l1,null,l2);
  }
  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    result = prime * result
        + ((characters == null) ? 0 : characters.hashCode());
    return result;
  }
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    PIPString other = (PIPString) obj;
    if (characters == null)
    {
      if (other.characters != null)
        return false;
    } else if (!characters.equals(other.characters))
      return false;
    return true;
  }
  
}