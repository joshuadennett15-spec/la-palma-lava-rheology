import pandas as pd
import sys
from pathlib import Path

file_path = Path(__file__).resolve().parent / 'LaPalmaData_Expanded_DATA_COMPREHENSION.xlsx'

try:
    xl = pd.ExcelFile(file_path)
    print('Sheet names:', xl.sheet_names)
    if len(xl.sheet_names) > 1:
        print('\n--- Second Sheet (' + xl.sheet_names[1] + ') ---')
        df = xl.parse(xl.sheet_names[1])
        print(df.to_string())
    else:
        print('Only one sheet found.')
except Exception as e:
    print('Error:', e)
