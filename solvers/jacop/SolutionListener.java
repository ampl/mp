import JaCoP.core.Var;
import JaCoP.search.SimpleSolutionListener;
import JaCoP.search.Search;
import JaCoP.search.SelectChoicePoint;

class SolutionListener<T extends Var> extends SimpleSolutionListener<T> {
  private long data;

  private static native void handleSolution(long data);

  SolutionListener(long data) { this.data = data; }
  
  public boolean executeAfterSolution(Search<T> search, SelectChoicePoint<T> select) {
    handleSolution(data);
    return super.executeAfterSolution(search, select);
  }
}
