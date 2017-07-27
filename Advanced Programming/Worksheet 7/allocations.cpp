#include <cstdlib>
class Pi {
  public:
    Pi(): v(new double){*v=3.1415;}
    ~Pi(){}
    double get() const {return *v;}
  private:
    double *v;
};

int main(){
  Pi a[3];
  Pi *b = new Pi[3];
  Pi *c = (Pi*) malloc(sizeof(Pi)*3);

  // TODO: remove potential memory leaks

  return 0;
}
