
import streamlit as st

from scripts.python.modelreader import ReadExplorerModel
from scripts.python.matcher import MatchSubmodel

# To work with local files in st 1.30.1, see
# https://discuss.streamlit.io/t/axioserror-request-failed-with-status-code-403/38112/13.
# The corresponding settings should not be used on a server.
uploader = st.sidebar.file_uploader("Model file (JSONL)")

# You can use a column just like st.sidebar:
srch = st.sidebar.text_input('Search pattern:')

fwd = st.sidebar.checkbox('Add descendants', disabled=True)
bwd = st.sidebar.checkbox('Add ancestors', disabled=True)

left_column, right_column = st.columns(2)


# Cache the reading function
@st.cache_data
def ReadModel(uploader):
  return ReadExplorerModel(uploader)

# Cache the matching function?
# @st.cache_data  Need cacheable Model.
def MatchSelection(m, srch, fwd, bwd):
  return MatchSubmodel(m, srch, fwd, bwd)

# Write dictionary of entries
def WriteDict(d):
  for k, v in d.items():
    if len(v):
      with st.expander("""### """ + k + ' (' + \
          str(v.count('\n')) + ')'):
            with st.container(height=200):
              st.code(v, language='ampl')

# Or even better, call Streamlit functions inside a "with" block:
if uploader is not None:
  model = ReadModel(uploader)
  subm1, subm2 = MatchSelection(model, srch, fwd, bwd)
  bytes1_data = subm1.GetData()
  bytes2_data = subm2.GetData()
  with left_column:
    st.write("""## NL model""")
    WriteDict(bytes1_data)
  with right_column:
    st.write("""## Flat model""")
    WriteDict(bytes2_data)
else:
  with left_column:
    st.write("No file selected.")
