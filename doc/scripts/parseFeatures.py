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


# Create the dataframe
df = create_dataframe(SOLVERS_PATH, feature_names, solver_names)

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
for v in backend_names_mapping.values():
    print(v)
print()
print(dataframe_to_rst(df))
