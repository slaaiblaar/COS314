import pandas as pd
import os

# Define the directories
dir1 = './results_static_mrate'
dir2 = './results_dynamic_mrate'

# Get the list of csv files in each directory
files1 = [f for f in os.listdir(dir1) if f.endswith('generations.csv')]
files2 = [f for f in os.listdir(dir2) if f.endswith('generations.csv')]

# Make sure the files are in the same order
files1.sort()
files2.sort()

# Loop through each file
for f1, f2 in zip(files1, files2):
    # Read the csv file from each directory
    df1 = pd.read_csv(os.path.join(dir1, f1))
    df2 = pd.read_csv(os.path.join(dir2, f2))

    # Drop the 'Generation' column from the second dataframe
    df2 = df2.drop(columns=['Generation'])

    # Rename the columns to include the directory name
    df1.columns = [f'Static {col}' if col != 'Generation' else col for col in df1.columns]
    df2.columns = [f'Dynamic {col}' for col in df2.columns]

    # Concatenate the dataframes horizontally
    result = pd.concat([df1, df2], axis=1)

    # Create a new Excel writer object
    writer = pd.ExcelWriter(f1[:-4] + '.xlsx')

    # Write the result to the Excel file
    result.to_excel(writer, sheet_name=f1[:-4], index=False)

    # Save the Excel file
    writer.close()
