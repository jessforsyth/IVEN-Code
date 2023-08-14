import pandas as pd
import numpy as np
from tkinter import filedialog	 #library to allow interactive file selection

print('\nExtract data from IVEN2output -')
print('Select the IVEN2out_.... file you would like to extract ')
filepath = filedialog.askopenfilename()

# load in data from excel file
data = pd.read_csv(filepath, header=[0])

# data is what is called a DataFrame, access any of the information by using the column header e.g.
data['IDs']   # gives the ID column
data['num_nbrs']  # gives the number of neighbours column

# then you can start to compare if you want