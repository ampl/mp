#ifndef BARONMPCOMMON_H
#define BARONMPCOMMON_H

#include "mp/backend-to-model-api.h"



#include <string>
#include <list>
#include <vector>
#include <variant>
#include <utility>
#include <memory>
#include <map>
#include <optional>

#include "mp/format.h"
#include "mp/posix.h" // for BufferedFile


#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace mp {



#ifdef _WIN32
  const std::vector<std::string> cplexlnames = { "cplex2212.dll", "cplex2211.dll","cplex2210.dll", "cplex2010.dll", "cplex12100.dll" };
  const std::vector<std::string> xpresslnames = { "xprs.dll" };
#elif defined(__APPLE__)
  const std::vector<std::string> cplexlnames = { "libcplex2211.dylib","libcplex2211.dylib", "libcplex2210.dylib", "libcplex2010.dylib","libcplex12100.dylib", "cplex.dylib" };
  const std::vector<std::string> xpresslnames = { "libxprs.dylib" };
#else
  const std::vector<std::string> cplexlnames = { "libcplex2212.so", "libcplex2211.so","libcplex2210.so", "libcplex2010.so","libcplex12100.so", "libcplex.so" };
  const std::vector<std::string> xpresslnames = { "libxprs.so.42", "libxprs.so.41", "libxprs.so.39.01", "libxprs.so.37.01", "libxprs.so.36.01.11","libxprs.so.36.01", "libxprs.so.34.01", "libxprs.so.33.01", "libxprs.so.32.01", "libxprs.so.31.01", "libxprs.so.29.01", "libxprs.so" };
#endif




  class Options {
  public:

    class OptionMetadata {
    public:
      using VT = std::variant<bool, int, double, std::string>;
      std::string key;
      std::string description;
      VT def;

      OptionMetadata(const std::string& key, const std::string& description,
        VT def)
        : key(key), description(description), def(def) {}
    };
    bool barstats = false;
    double deltaa = INFINITY;
    int deltaterm = 0;
    double deltar = 1;
    double deltat = -100.0;
    double epsa = 1e-6;
    double epsr = 1e-9;

    int firstfeas = 0;
    int firstloc = 0;

    int iismethod = 0;
    int iisint = 0;
    int iisorder = 0;
    bool keepsol = false;

    // Sorted by ValidateLPSolver()
    int lpsol;
    std::string cplexlibname = "";
    std::string xpresslibname = "";
    // end

    std::string lpsolver = "cbc";
    bool lsolmsg = false;
    std::string lsolver;

    int maxiter = -1;
    double maxtime = 1000.0;

    int nlpsolver = 0;
    int numsol = 1;

    bool objbound = false;
    int objno = 1;

    int optfile = 0;
    std::string output_filename;

    int outlev = 0;

    double prfreq = 1000000;
    int prloc = 0;
    std::string problem = "amplproblem";
    double prtime = 30.0;
    std::string scratch = "";
    int seed = 19631963;
    std::string sumfile = "";
    int threads = 1;
    std::string trace = "";
    // Not serialized
    bool overwrite = false;

    std::vector<std::string> inlineparams;

    void initializeLPSolver();
    // Map from field names to OptionMetadata
    const std::map<std::string, OptionMetadata> option_metadata = {

        {"alg:deltaa deltaa", {"DeltaA", "Adjustment parameter for delta values.", INFINITY}},
        {"alg:deltaterm deltaterm", {"DeltaTerm", "Delta termination criterion.", 0}},
        {"alg:deltar deltar", {"DeltaR", "Relative delta adjustment.", 1.0}},
        {"alg:deltat deltat", {"DeltaT", "Time-based delta adjustment.", -100.0}},
        {"alg:epsa epsa", {"EpsA", "EpsA convergence tolerance (default 1e-6). "
            "BARON stops if the current function value f satisfies "
            "abs(f - L) <= epsa, where L is the currently best available bound on f.", 1e-6}},
        {"alg:epsr epsr", {"EpsR", "BARON's EpsR convergence tolerance (default 1e-9). "
            "BARON stops if the current function value f satisfies "
            "abs(f - L) <= abs(L*epsr), where L is the currently best available bound on f.", 1e-9}},
        {"alg:firstfeas firstfeas", {"FirstFeas", "If set to 1, BARON will terminate once it finds "
            "numsol feasible solutions, irrespective of solution quality. Default is 0, "
            "meaning that BARON will search for the *best* numsol feasible solutions.", 0}},
        {"alg:firstloc firstloc", {"FirstLoc", "If set to 1, BARON will terminate once it finds "
            "a local optimum, irrespective of solution quality. Default is 0, "
            "meaning that BARON will search for the best numsol feasible solutions.", 0}},
       // TODO
       /* {"alg:iisint iisint", {"IISint", "Whether to include integer variables in an IIS:\n\n.. value-table::\n",0}},
        {"alg:iismethod iismethod", {"CompIIS", "Which method to use when finding an IIS (irreducible infeasible "
              "set of constraints, including variable bounds):\n\n.. value-table::\n", 0}},
        {"alg:iisorder iisorder", {"IISorder","How to order constraints when seeking an IIS:\n\n.. value-table::\n",0}},*/



      // Not serialized - at least not trivially
      {"tech:barstats barstats", {"", "Report detailed BARON statistics.", false}},
      {"tech:keepsol keepsol", {"", "Keep BARON's solution files.", false}},
      {"tech:lpsolver lpsolver", {"", "Choice of LP solver, which matters mainly when there are integer "
                                      "variables:  one of cbc (default), cplex, or xpress.  The last two "
                                      "must be suitably licensed to be used.", ""}},
      {"tech:cplexlibname cplexlibname", {"CplexLibName", "If used, specifies the path to cplex  libraries", ""}},
      {"tech:xpresslibname xpresslibname", {"XpressLibName", "If used, specifies the path to xpress libraries", ""}},
      {"tech:lsolver lsolver", {"", "Local nonlinear solver that BARON should call.", ""}},
      {"tech:lsolmsg lsolmsg", {"", "Show solution messages for lsolver.", ""}},
      {"tech:nlpsolver nlpsol nlpsolver", {"", "Local nonlinear solvers BARON is allowed to use: sum (mod 16) of\n"
          "| 1 - IPOPT (builtin)\n"
          "| 2 - FilterSD (builtin)\n"
          "| 4 - FilterSQP (builtin)\n"
          "| 8 - lsolver (if lsolver=... is specified)\n"
          "Default 0 - allow all.", 0}},

        {"tech:numsol numsol", {"NumSol", "Number of near optimal solutions to find. "
                                  "Default = 1; values > 1 imply keepsol and cause suffix .numsol "
                                  "on the objective and problem to be returned.", 1}},
        {"lim:iter iterlimit maxiter", {"MaxIter", "Maximum number of iterations.", -1}},
        {"lim:time timelim timelimit maxtime", {"MaxTime", "Maximum time in seconds (default 1000).", 1000.0}},


        {"tech:objbound objbound", {"", "Return suffixes .obj_lb and .obj_ub on the problem and objective with "
                                        "Baron's final lower and upper bounds on the objective value", false}},
        {"tech:objno objno", {"", "Objective number (default: 1).", 1}},
        {"tech:optfile optfile", {"", "Name of BARON option file (not required).  If given, the file should "
                                "contain name-value pairs, one per line, with the name and value "
                                "separated by a blank, a colon, or an equal sign, possibly surrounded "
                                "by white space.  The names and possible values are summarized in "
                                "section 7 of the BARON user manual (baron_manual.pdf).  Empty lines "
                                "and lines that start with # are ignored.", ""}},
        {"tech:outlev outlev", {"", "Output verbosity level.", 0}},
        {"tech:overwrite overwrite", {"", "If set, overwrite the files in the scratch directory.", 0}},

        {"tech:prfreq prfreq", {"PrFreq", "Report progress every prfreq nodes (default 1e6).", 1e6}},
        {"tech:prloc prloc", {"LocRes", "Whether to report local searches:\n\n.. value-table::\n", 0}},
        {"tech:problem problem", {"ProName", "Problem name printed in logfile.", "amplproblem"}},
        {"tech:prtime prtime", {"PrTimeFreq", "Report progress every prtime seconds  (default 30).", 30.0}},
        {"tech:scratch scratch", {"scratch", "Directory for temporary files; if set to 'local' it uses a subdirectory of the NL file location with the NL file name", std::string("")}},
        {"tech:seed seed", {"seed", "Initial seed for random number generation, must be a positive integer (default 19631963).", 19631963}},
        {"tech:sumfile sumfile", {"SumName", "Name of summary file.", std::string("")}},
        {"tech:threads threads", {"threads", "Maximum number of threads to use (default 1 when integer variables are present)", 1}},
        {"tech:trace trace", {"trace", "Name of BARON trace file.", std::string("")}}
    };

    void setProblemName(std::string_view stub);
    const char* description(std::string_view name) {
      return option_metadata.at(name.data()).description.c_str();
    }

    void serializeMember(std::string_view name, OptionMetadata::VT v, fmt::MemoryWriter& m, bool alwaysSerialize=false) const {
      
      auto md = option_metadata.at(name.data());
      if (alwaysSerialize || (md.def != v)) {
        m << md.key << ": ";
        std::visit([&m](auto&& arg) {
          if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, std::string>) {
            m << "\"" << arg << "\"";  // Surround strings with quotes
          }
          else {
            m << arg;
          }
          }, v);
        m << ";\n";
      }
    }

#define ADDSER(var, n)\
      serializeMember(n, var, m)

#define ADDALGSER(n) ADDSER(n, "alg:"#n" "#n)
#define ADDTECHSER(n) ADDSER(n, "tech:"#n" "#n)


    // Serialization method using fmt::MemoryWriter
    std::string serialize() const {
      fmt::MemoryWriter m;
      ADDALGSER(deltaa);
      ADDALGSER(deltaterm);
      ADDALGSER(deltar);
      ADDALGSER(deltat);
      ADDALGSER(epsa);
      ADDALGSER(epsr);
      ADDALGSER(firstfeas);
      ADDALGSER(firstloc);
     /* ADDALGSER(iisint);
      ADDALGSER(iismethod);
      ADDALGSER(iisorder);*/
      ADDSER(maxiter, "lim:iter iterlimit maxiter");
      ADDSER(maxtime, "lim:time timelim timelimit maxtime");
      ADDTECHSER(numsol);
      //ADDTECHSER(optfile);
      ADDTECHSER(prfreq);
      serializeMember("tech:prloc prloc", prloc, m, true);
      ADDTECHSER(problem);
      ADDTECHSER(prtime);
      ADDTECHSER(seed);
      ADDTECHSER(sumfile);
      ADDTECHSER(threads);
      ADDTECHSER(trace);
      m << fmt::format("LPSol: {};\n", lpsol);
      if( (lpsol==3) && (!cplexlibname.empty()))
          ADDTECHSER(cplexlibname);
      // ADDTECHSER(xpresslibname); NOT supported anymore?
      appendInlineParams(m);
      m << fmt::format("PrLevel: {};\n", outlev);
      return m.str();
    }
    void appendInlineParams(fmt::MemoryWriter& m) const;
  };

  




  /// Information shared by both
/// `BaronmpBackend` and `BaronmpModelAPI`
struct BaronmpCommonInfo {
  struct BaronProblemInfo {
    std::vector<std::string> varNames;
    std::vector<std::string> conNames;
    std::map<std::string, int> varMap;
    std::vector<int> baronToAMPLIndices;
    int nVarsInteger = 0, nVarsBinary = 0, nVarsContinuous = 0;
    int obj_sense;
  };
  std::shared_ptr <Options> baronOptions_;
  Options& baronOptions() {
    return *baronOptions_;
  }
  
  int nVarsInteger()    const { return lp()->nVarsInteger;} 
  int nVarsBinary()     const{ return lp()->nVarsBinary;} 
  int nVarsContinuous() const {return lp()->nVarsContinuous; } 
  void nVarsInteger(int n)    {  lp()->nVarsInteger=n; }
  void nVarsBinary(int n)     {  lp()->nVarsBinary=n;}
  void nVarsContinuous(int n) { lp()->nVarsContinuous=n; }
  
  // Directory where to write baron files
  std::string baronDir;
  // NL file path - shouldn't really need it but let's see TODO
  std::string nlFilePath;
  // Directory where the execution has initiated
  std::string initialDir;
  std::string filePathBar_, filePathDic_, filePathAMPL_;
  std::string filePathBar() {
    return filePathBar_;
  }
  std::string filePathDic() {
    return filePathDic_;
  }
  std::string filePathAMPL() {
    return filePathAMPL_;
  }
  std::shared_ptr<fmt::BufferedFile> FILE_BAR, FILE_DIC, FILE_AMPL;

  BaronProblemInfo* lp() const { return lp_; }
  void set_lp(BaronProblemInfo* lp) { lp_ = lp; }

  BaronmpCommonInfo() {
    baronOptions_ = std::make_shared<Options>(Options());
  }
private:
  BaronProblemInfo*      lp_ = NULL;

};




/// Common API for Baronmp classes
class BaronmpCommon :
    public Backend2ModelAPIConnector<BaronmpCommonInfo> {
public:

#ifndef WIN32
  static volatile pid_t pid;
  #else 
  static volatile DWORD pid;
#endif
  // Define version as follows
  const int v_day =10, v_month = 12, v_year = 25;
  int currentObj = 0;
  static constexpr double Infinity() { return INFINITY;  }
  static constexpr double MinusInfinity() { return -INFINITY; }
  const std::string FILENAME_BAR = "amplmodel";
  const std::string FILENAME_DIC = "dictionary";
  const std::string FILENAME_AMPL = ".amplparams";
  const std::string BARON_TIM = "tim.lst";
  const std::string BARON_RES = "res.lst";



  std::string appendToDir(const std::string& dir, const std::string& file);

  void writeBaronOptions();
  void writeBaron(const std::string &s) {
    std::fwrite(s.c_str(), s.size(), 1, FILE_BAR->get());
  }
  void writeVars(fmt::MemoryWriter &w, const std::string& s) {
    w << s;
  }
  void writeVars(fmt::MemoryWriter& w, fmt::CStringRef format, const fmt::ArgList& args) {
    w << fmt::format(format, args);
  }
  /// Formats a string and prints it to stdout or, if an output handler
/// is registered, sends it to the output handler.
  void writeBaron(fmt::CStringRef format, const fmt::ArgList& args) {
    FILE_BAR->print(format, args);
  }
  /// Variadic overload of Print()
  FMT_VARIADIC(void, writeBaron, fmt::CStringRef)

  int runBaron( const std::string& arg, double timelimit);
  std::string make_cmdline(const std::vector<std::string>& args);
  int run(const std::vector<std::string>& args, double timelimit);
  void initDirectories(const std::string& stub, const std::string& scratch, bool overwrite);

  void initBaronFile();
  void deinitBaronFile(bool changedir=true);
  int recrmdir(const std::string& dname);
  void changeDirectory(const std::string& path);
  
};

class BaronGlobals {
  public:
    std::string lsolver = "[builtin]";
    int lsolmsg;           // Expect messages for subsolver
    int nvars;             
    std::string nlfile;    // from getNames
    int v_keepsol;         // from option
    std::string initialDir;
    std::string baronDir;
    std::string verbuf;
    std::string lpsol_dll;

  // Write the globals data to a file
  bool serialize(const std::string& filename) ;
  // Read the globals data from a file
  static std::optional<BaronGlobals> deserialize(const std::string& filename) ;
};


class TimFileData {
public:
  double printlb;
  double printub;
  int barstatus;
  int modelstatus;
  int badmissingbounds;
  int itera;
  int nodeopt;
  int nodememmax;
  double totaltime;
  int Error_code;

  // Default constructor
  TimFileData()
    : printlb(0), printub(0), barstatus(0), modelstatus(0),
    badmissingbounds(0), itera(0), nodeopt(0), nodememmax(0),
    totaltime(0.0), Error_code(0) {}

  // Static factory method to create an instance from file
  static TimFileData fromFile(const std::string& filePath); 
private:
 
};



class ResFileData {


  // Member variables to store the parsed data
  std::string baronVersion = "BARON version ";
  std::string baronLP = " LP solver:";
  std::string baronNLP = " NLP solver:";
  std::string baronNumsol = "numsol         =";
  std::string baron_version = "BARON version ";
  std::string baron_terminated_solution = "The above solution has an objective value of:";
  std::string baron_terminated_start = "variable";
  std::string baron_objective = "Objective value is:";
  std::string baron_dual = "Corresponding dual";

  std::vector<double> objective_values;
  std::vector<std::vector<double>> primal_solutions;
  std::vector<std::vector<double>>  dual_solutions;
  bool baron_terminated;
  bool baron_solution_read;
  int pos_betterfound;
  int numsol = 0;
  int Error_code;
  
public:

  int NumSol() { return numsol; }
  const std::vector<double>& Objectives() { return objective_values; }


  double ObjectiveValue() const {
    
    if(objective_values.size() == 0)
    return 0;
    return objective_values[objective_values.size() - 1];
  }

  const std::vector<std::vector<double>>& PrimalSolutions() {
    return primal_solutions;
  }
  const std::vector<std::vector<double>>& DualSolutions() {
    return dual_solutions;
  }
  const std::vector<double>& PrimalSolution() {
    if (primal_solutions.size() > 0)
      return primal_solutions[primal_solutions.size() - 1];
    else
    {
      primal_solutions.push_back(std::vector<double>());
      return primal_solutions[0];
    }
  }
  const std::vector<double>& DualSolution() {
    if (dual_solutions.size() > 0)
        return dual_solutions[dual_solutions.size() - 1];
    else
    {
      dual_solutions.push_back(std::vector<double>());
      return dual_solutions[0];
    }
  }
  // Default constructor
  ResFileData() : baron_terminated(false), baron_solution_read(false),
    pos_betterfound(-1), numsol(0), Error_code(0) {}

  // Static factory method to create an instance from file
  static ResFileData fromFile(const std::string& filePath, int nodeopt);

private:
  // Helper function to parse baron-specific data in a line
  void parseBaronInfo(const std::string& line);
  bool readSolution(std::ifstream& barRes);
};


class Utils {
public:
  static std::string getLibPath();
  /**
  * Baron can output scientific notation with D instead of E
  */
  static std::optional<double> readBaronValues(std::string value)
  {
    // Replace first occurence
    auto pos = value.find_first_of("D");
    if (pos != std::string::npos) {
      value[pos] = 'E';
    }
    try {
      return std::stod(value);
    }
    catch (const std::exception&) {
      return std::nullopt;
    }
  }
  static std::string findFile(const std::string& filename, bool isDLL = false);
private:
  static std::string getMyPath(bool getDir, const std::string& appendPath = "");
  static std::string getExecutablePath();
 
 
};

} // namespace mp

#endif // BARONMPCOMMON_H
