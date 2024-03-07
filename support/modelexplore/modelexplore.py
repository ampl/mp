
import streamlit as st

# To work in st 1.30.1, see
# https://discuss.streamlit.io/t/axioserror-request-failed-with-status-code-403/38112/13
input_file = st.file_uploader("Model file (JSONL)")

left_column, right_column = st.columns(2)

# You can use a column just like st.sidebar:
srch = left_column.text_input('Search pattern:')

# Or even better, call Streamlit functions inside a "with" block:
fwd = right_column.checkbox('Add descendants', disabled=True)
bwd = right_column.checkbox('Add ancestors', disabled=True)

if input_file is not None:
  with left_column:
    bytes_data = input_file.read()
    st.write("NL model")
    st.write(bytes_data)
  with right_column:
    st.write("Flat model")
else:
  with left_column:
    st.write("No file selected.")
