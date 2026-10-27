"""Feature support of the MP-based drivers, read from their ALLOW_STD_FEATURE
declarations.

Without arguments, prints the support table for the MP features guide
(doc/source/features-guide.rst, section "Support by solvers").
With --devampl DIR, also writes one Markdown fragment per solver,
DIR/<solver>/features.inc, which the solver pages on dev.ampl.com include.
"""
import argparse
import re
import os
import pandas as pd

default_features = {
   "MULTIOBJ" : 2,
   "REPORT_TIMES" : 1
}
backend_names_mapping = {}

LOCAL_PATH =  os.path.dirname(os.path.abspath(__file__))
SOLVERS_PATH = os.path.join(LOCAL_PATH, "..", "..", "solvers")
# Path to the feature names CSV file
feature_names = os.path.join(LOCAL_PATH, "features_names.txt")
solver_names = os.path.join(LOCAL_PATH,"solver_names.txt")

# Function to parse ALLOW_STD_FEATURE macros in a C++ file
def parse_macros(file_path):
    with open(file_path, "r") as file:
        content = file.read()
    
    # Pattern to match ALLOW_STD_FEATURE macros
    macro_pattern =re.compile(r"ALLOW_STD_FEATURE\s*\(\s*(\w+)\s*,\s*(true|false)\s*\)")
    macros = macro_pattern.findall(content)
    return {feature: value == "true" for feature, value in macros}

# Function to find all _backend.h files recursively in a directory
def find_backend_files(directory):
    backend_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith("backend.h"):
                backend_files.append(os.path.join(root, file))
    return backend_files

# Function to read the features_names.csv file
def read_feature_names(file_path):
    feature_names = {}
    with open(file_path, "r") as file:
        for line in file:
            macro_name, friendly_name = line.strip().split(",")
            feature_names[macro_name] = friendly_name
    return feature_names

def extract_ref(rstlink: str):
    # .. _Copt: https://dev.ampl.com/solvers/copt/
    macro_pattern =re.compile(r".. _(\w+):")
    macros = macro_pattern.findall(rstlink)
    return macros[0]+"_"

# Main function to create the dataframe
def create_dataframe(directory, feature_names_file, backend_mapping_file):
    # Find all backend files
    backend_files = find_backend_files(directory)

    # Read feature names from the CSV file
    feature_names = read_feature_names(feature_names_file)

    # Read the order of features from the CSV file
    with open(feature_names_file, 'r') as file:
        feature_order = [line.strip().split(',')[0] for line in file]

    # Read backend mapping from the CSV file
    global backend_names_mapping
    with open(backend_mapping_file, 'r') as file:
        for line in file:
            backend_name, friendly_name = line.strip().split(',')
            backend_names_mapping[backend_name] = friendly_name

    # Filter backend names based on the mapping file
    backendfile_names = [os.path.basename(file_path) for file_path in backend_files] 
    filtered_backends = [backend for backend in backend_names_mapping.keys() if backend in backendfile_names]

    # Dictionary to store macros for each backend file
    all_macros = {}
    for file_path in backend_files:
        if  os.path.basename(file_path) in filtered_backends:
            all_macros[os.path.basename(file_path)] = parse_macros(file_path)

    # Create a dataframe
    data = []
    for feature in feature_order:
        row = [feature_names.get(feature, feature)]  # Use friendly name if available
        for backend_name in filtered_backends:
            value = all_macros.get(backend_name, {}).get(feature, default_features.get(feature,False))
            # MP emulates multiple objectives for every driver without native support
            if feature == "MULTIOBJ" and value is False:
                value = 2
            row.append('|y|' if value==1 else '|n|' if value ==0 else '|e|')
        data.append(row)

    columns = ['Feature'] + [extract_ref(backend_names_mapping[backend]) for backend in filtered_backends]
    df = pd.DataFrame(data, columns=columns)

    return df


def dataframe_to_rst(df):
    # Get column names and maximum lengths
    columns = list(df.columns)
    max_lengths = [max(df[column].apply(len)) for column in columns]
    length_headers = [len(h) for h in columns]
    max_lengths = [max(c, h) for c,h in zip(max_lengths,length_headers)]
    separator = "+" + "+".join(["-" * (length + 2) for length in max_lengths]) + "+\n"
    # Create RST table header
    header = separator
    header += "|" + "|".join([f" {column.ljust(length)} " for column, length in zip(columns, max_lengths)]) + "|\n"
    header += "+" + "+".join(["=" * (length + 2) for length in max_lengths]) + "+\n"
    
    # Create RST table rows
    rows = ""
    for index, row in df.iterrows():
        values = [value for value in row]
        rows += "|" + "|".join([f" {value.ljust(length)} " for value, length in zip(values, max_lengths)]) + "|\n"
        rows+=separator
    
    # Combine header and rows
    rst_table = header + rows
    
    return rst_table
MP_FEATURES = "https://mp.ampl.com/features-guide.html"

# Feature -> (title on dev.ampl.com, main option, anchor in the MP features guide)
DEVAMPL_FEATURES = {
    "WRITE_PROBLEM": ("Model export", "tech:writemodel", "model-export"),
    "WRITE_SOLUTION": ("Solution export", "tech:writesolution", None),
    "FEAS_RELAX": ("Feasibility relaxation", "alg:feasrelax", "feasibility-relaxation"),
    "IIS": ("Irreducible inconsistent subset (IIS)", "iis:find",
            "irreducible-inconsistent-subset-iis"),
    "BASIS": ("Basis input and output", "alg:basis", "input-and-output-basis"),
    "KAPPA": ("Condition number (kappa)", "alg:kappa", "kappa"),
    "MULTISOL": ("Multiple solutions", "sol:stub", "multiple-solutions"),
    "MULTIOBJ": ("Multiple objectives", "obj:multi", "multiple-objectives"),
    "SENSITIVITY_ANALYSIS": ("Sensitivity analysis", "alg:sens", "sensitivity-analysis"),
    "RAYS": ("Unbounded rays and Farkas certificates", "alg:rays", "unbounded-rays"),
    "WARMSTART": ("Warm start", "alg:start", "warm-start"),
    "FIX_MODEL": ("Fixed model (duals and basis for a MIP)", "mip:basis",
                  "fixed-model-return-basis-for-mip"),
    "LAZY_USER_CUTS": ("Lazy constraints and user cuts", "mip:lazy",
                       "lazy-constraints-and-user-cuts"),
    "RETURN_MIP_GAP": ("MIP gap", "mip:return_gap", "return-mip-gap"),
    "RETURN_BEST_DUAL_BOUND": ("Best dual bound", "mip:bestbound", "return-best-dual-bound"),
    "VAR_PRIORITIES": ("Branching priorities", "mip:priorities", "variable-priorities"),
    "PLATEAU_STOP": ("Stop on a plateau", "mip:plateau:time",
                     "plateau-stop-stall-based-mip-termination"),
    "PLATEAU_STOP_BOUND": ("Stop on a plateau of the MIP gap", "mip:plateau:relgaptol",
                           "plateau-stop-including-mip-gap-progress"),
}
# Text for features that have no section in the MP features guide
DEVAMPL_NO_ANCHOR = {
    "WRITE_SOLUTION": "Writes the solution in the solver's native format, "
                      "e.g., `writesol=model.sol`.",
}
MIP_ONLY = {"FIX_MODEL", "LAZY_USER_CUTS", "RETURN_MIP_GAP", "RETURN_BEST_DUAL_BOUND",
            "VAR_PRIORITIES", "PLATEAU_STOP", "PLATEAU_STOP_BOUND"}

# Solver particulars shown under a feature: (solver page, feature) -> Markdown.
# Option names must exist in the driver's -= table.
DEVAMPL_NOTES = {
    ("gurobi", "IIS"): "Method: `iis:method`; IIS even for a feasible model "
                       "part: `iis:force`.",
    ("gurobi", "KAPPA"): "Exact computation: `alg:kappa_exact`.",
    ("gurobi", "FEAS_RELAX"): "Big-M for relaxing indicator constraints: "
                              "`alg:feasrelaxbigm`.",
    ("gurobi", "MULTIOBJ"): "Native hierarchical and blended objectives; "
                            "`obj:multi:method`, `obj:multi:inheritparams`.",
    ("gurobi", "FIX_MODEL"): "Algorithm for the fixed model: `mip:fixedmethod`.",
    ("xpress", "IIS"): "Number of IISs: `iis:max`; scope: `iis:ops`.",
    ("cplex", "FIX_MODEL"): "Algorithm for the fixed model: `mip:basis:method`.",
    ("copt", "IIS"): "Method: `iis:method`.",
    ("highs", "IIS"): "Strategy: `iis:strategy`.",
    ("cpsat", "WARMSTART"): "Initial values are passed as CP-SAT solution hints.",
    ("mp2nl", "*"): "The features below are handled by MP2NL itself; others "
                    "depend on the subsolver and its options.",
}


def devampl_page(url):
    m = re.search(r"dev\.ampl\.com/solvers/(\w+)/", url)
    return m.group(1) if m else None


def devampl_fragment(page, support):
    """Markdown for one solver: a dropdown per supported feature."""
    lines = ["<!-- Generated by mp/doc/scripts/parseFeatures.py from the driver's",
             "     feature declarations. Do not edit: rerun the script. -->", ""]
    if (page, "*") in DEVAMPL_NOTES:
        lines += [DEVAMPL_NOTES[(page, "*")], ""]
    for heading, mip in (("For all models", False), ("For MIP models", True)):
        items = [f for f in DEVAMPL_FEATURES
                 if (f in MIP_ONLY) == mip and support.get(f, 0)]
        if not items:
            continue
        lines += [f"**{heading}**", ""]
        for f in items:
            title, option, anchor = DEVAMPL_FEATURES[f]
            body = f"Option `{option}`."
            if anchor:
                body += f" [Details]({MP_FEATURES}#{anchor})"
            else:
                body += " " + DEVAMPL_NO_ANCHOR[f]
            if support[f] == 2:
                body = "Emulated by MP. " + body
            lines += [f":::{{dropdown}} {title}", body]
            if (page, f) in DEVAMPL_NOTES:
                lines += ["", DEVAMPL_NOTES[(page, f)]]
            lines += [":::", ""]
    return "\n".join(lines)


def write_devampl(directory, solvers_dir):
    for backend, rstlink in backend_names_mapping.items():
        page = devampl_page(rstlink)
        files = [f for f in find_backend_files(solvers_dir)
                 if os.path.basename(f) == backend]
        if not page or not files:
            continue
        macros = parse_macros(files[0])
        support = {}
        for f in DEVAMPL_FEATURES:
            v = macros.get(f, default_features.get(f, False))
            if f == "MULTIOBJ" and v is False:
                v = 2
            support[f] = int(v)
        os.makedirs(os.path.join(directory, page), exist_ok=True)
        path = os.path.join(directory, page, "features.inc")
        with open(path, "w") as out:
            out.write(devampl_fragment(page, support))
        print("wrote", path)


parser = argparse.ArgumentParser(
    description="Feature support of the MP-based drivers, from their "
                "ALLOW_STD_FEATURE declarations.",
    epilog="Examples:\n"
           "  python parseFeatures.py   # print the table for features-guide.rst\n"
           "  python parseFeatures.py --devampl ../../../ampl.github.io/docs/source/solvers",
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--devampl", metavar="DIR",
                    help="also write DIR/<solver>/features.inc for each solver whose "
                         "link in solver_names.txt is a dev.ampl.com solver page; "
                         "DIR is docs/source/solvers of the ampl.github.io repository. "
                         "Example: --devampl ~/ampl.github.io/docs/source/solvers")
args = parser.parse_args()

df = create_dataframe(SOLVERS_PATH, feature_names, solver_names)
for v in backend_names_mapping.values():
    print(v)
print()
print(dataframe_to_rst(df))
if args.devampl:
    write_devampl(args.devampl, SOLVERS_PATH)
