import pandas as pd
import requests
import json
from pandas import json_normalize
import time

# Constants
chunk_size = 200
sleep_time = 60  # in seconds
input_file_path = './SRTR_input_for_eplet_API_UPDATED.csv'
output_dir = './SRTR_output_API/'

# List to store incomplete chunks
incomplete_chunks = []
failed_chunks = []

# Read in api key file, saved locally and not checked in
EpReg_api_key_file = open ("HLA_EpReg_api_key.key")
EpReg_api_key = ""
for line in EpReg_api_key_file:
    EpReg_api_key = line

# Function to process each chunk and log errors
def process_chunk(chunk, chunk_number):
    eplet_list = []
    # 1st Error handling for chunk processing 
    try:
      for index, row in chunk.iterrows():
          # 2nd Error handling for requests and JSON decoding
          try:
              # Extract donor and recipient information from the row
              PX_ID = row['PX_ID']
              DONOR_A_1 = row["DONOR_A_1"]
              DONOR_A_2 = row["DONOR_A_2"]
              DONOR_B_1 = row["DONOR_B_1"]
              DONOR_B_2 = row["DONOR_B_2"]
              DONOR_C_1 = row["DONOR_C_1"]
              DONOR_C_2 = row["DONOR_C_2"]
              DONOR_DRB1_1 = row["DONOR_DRB1_1"]
              DONOR_DRB1_2 = row["DONOR_DRB1_2"]
              DONOR_DRW1_1 = row["DONOR_DRB345_1"]
              DONOR_DRW1_2 = row["DONOR_DRB345_2"]
              DONOR_DQA1_1 = row["DONOR_DQA1_1"]
              DONOR_DQA1_2 = row["DONOR_DQA1_2"]
              DONOR_DQB1_1 = row["DONOR_DQB1_1"]
              DONOR_DQB1_2 = row["DONOR_DQB1_2"]
              RECIP_A_1 = row["RECIP_A_1"]
              RECIP_A_2 = row["RECIP_A_2"]
              RECIP_B_1 = row["RECIP_B_1"]
              RECIP_B_2 = row["RECIP_B_2"]
              RECIP_C_1 = row["RECIP_C_1"]
              RECIP_C_2 = row["RECIP_C_2"]
              RECIP_DRB1_1 = row["RECIP_DRB1_1"]
              RECIP_DRB1_2 = row["RECIP_DRB1_2"]
              RECIP_DRW1_1 = row["RECIP_DRB345_1"]
              RECIP_DRW1_2 = row["RECIP_DRB345_2"]
              RECIP_DQA1_1 = row["RECIP_DQA1_1"]
              RECIP_DQA1_2 = row["RECIP_DQA1_2"]
              RECIP_DQB1_1 = row["RECIP_DQB1_1"]
              RECIP_DQB1_2 = row["RECIP_DQB1_2"]
              #print(RECIP_DQB1_1)

              immunizer_string = DONOR_A_1+",,,"+ DONOR_A_2+",,,"+  DONOR_B_1+",,,"+ DONOR_B_2+",,,"+ DONOR_C_1+",,,"+ DONOR_C_2+",,,"+ DONOR_DRB1_1+",,,"+ DONOR_DRB1_2+",,,"+ DONOR_DRW1_1+",,,"+ DONOR_DRW1_2+",,,"+ DONOR_DQA1_1+",,,"+DONOR_DQA1_2+",,,"+DONOR_DQB1_1+",,,"+DONOR_DQB1_2

              #print (immunizer_string)

              patient_string = RECIP_A_1+",,,"+ RECIP_A_2+",,,"+ RECIP_B_1+",,,"+ RECIP_B_2+",,,"+ RECIP_C_1+",,,"+ RECIP_C_2+",,,"+ RECIP_DRB1_1+",,,"+ RECIP_DRB1_2+",,,"+ RECIP_DRW1_1+",,,"+ RECIP_DRW1_2+",,,"+ RECIP_DQA1_1+",,,"+ RECIP_DQA1_2+",,,"+ RECIP_DQB1_1+",,,"+ RECIP_DQB1_2

              #print (patient_string)

              request_path = f'https://api.epregistry.com.br/eplet_mismatches?from={EpReg_api_key}&immunizer_alleles={immunizer_string}&patient_alleles={patient_string}'

              # Call the API
              r = requests.get(request_path, headers={'Accept': 'application/json'}, timeout=30) # response object (contains JSON-formatted content)
              # read timeout increased to 30s
              r.raise_for_status()  # Check if the request was successful # Raise HTTPError for bad responses
              output = r.json() # converts json-formatted response object into a Python dictionary
              # Raises requests.exceptions.JSONDecodeError if response is not valid JSON

              #print ("API output: ")
              #print (output)

              # Flatten JSON-formatted nested data into row
              output_df = json_normalize(output, sep='_')

              replacements = {
                  DONOR_A_1: 'A_allele_1',
                  DONOR_A_2: 'A_allele_2',
                  DONOR_B_1: 'B_allele_1',
                  DONOR_B_2: 'B_allele_2',
                  DONOR_C_1: 'C_allele_1',
                  DONOR_C_2: 'C_allele_2',
                  DONOR_DRB1_1: 'DRB1_allele_1',
                  DONOR_DRB1_2: 'DRB1_allele_2',
                  DONOR_DRW1_1: 'DRB345_allele_1',
                  DONOR_DRW1_2: 'DRB345_allele_2',
                  DONOR_DQA1_1: 'DQA1_allele_1',
                  DONOR_DQA1_2: 'DQA1_allele_2',
                  DONOR_DQB1_1: 'DQB1_allele_1',
                  DONOR_DQB1_2: 'DQB1_allele_2'
              }

              # Create the new column names list
              col_names = []

              for col in output_df.columns:
                for old, new in replacements.items():
                  col = col.replace(old, new)
                col_names.append(col)

              # Reassign the new column names to the DataFrame
              output_df.columns = col_names

              # Include identifier
              output_df.loc[0, 'PX_ID'] = PX_ID

              # Append each dataframe output to the empty list
              eplet_list.append(output_df)

          except requests.exceptions.Timeout as e:
              error_message = f"Error processing row {index} with donor-recipient pair {row['PX_ID']} in chunk {chunk_number}: {e}"
              print(error_message)
              log_error(error_message)
              # EDIT: Implement a retry mechanism?
              incomplete_chunks.append(chunk_number)

          except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
              error_message = f"Error processing row {index} with donor-recipient pair {row['PX_ID']} in chunk {chunk_number}: {e}"
              print(error_message)
              log_error(error_message)
              incomplete_chunks.append(chunk_number)

              # Append a blank row with just the 'PX_ID' filled
              #blank_row = pd.DataFrame({col: [None] for col in output_df.columns})
              #blank_row['PX_ID'] = PX_ID
              #eplet_list.append(blank_row)
              #continue # Move to the next row within this chunk

      if eplet_list:
          eplet_dataframe = pd.concat(eplet_list, ignore_index=True)
          output_file_path = f"{output_dir}eplet_output_chunk_{chunk_number}.csv"
          eplet_dataframe.to_csv(output_file_path, index=False)
          print(f"Saved chunk {chunk_number} to {output_file_path}")

      return True # Chunk processed successfully

    except Exception as e:
        error_message = f"Error processing chunk {chunk_number}: {str(e)}"
        print(error_message)
        log_error(error_message)
        failed_chunks.append(chunk_number) # Mark the entire chunk as failed

        return False  # Chunk failed to process

    #return eplet_dataframe

# Main loop to process the dataset in chunks
#chunk_number = 0
#for chunk in pd.read_csv(input_file_path, chunksize=chunk_size):
#    chunk_number += 1
#    print(f"Processing chunk {chunk_number}...")
#    if process_chunk(chunk, chunk_number):
#      print(f"Processed chunk {chunk_number} successfully, sleeping for {sleep_time} seconds...")
#      time.sleep(sleep_time)
#    else:
#      print(f"First pass: failed to process chunk {chunk_number}, logging for further review and sleeping for {sleep_time} seconds...")
#      time.sleep(sleep_time)

# EDIT: Indicate which run # to log errors
# Function to log errors to a file
def log_error(message):
    with open("error_log_run3.txt", "a") as log_file:
        log_file.write(f"{message}\n")

# To rerun specific chunks, you can read the incomplete_chunks.json and process those chunks again
def rerun_incomplete_chunks(incomplete_chunks_file):
    with open(incomplete_chunks_file, "r") as f:
        incomplete_chunks = json.load(f)

    for chunk_number in incomplete_chunks:
        skip_rows = (((chunk_number - 1) * chunk_size) + 1) # ADJUSTED from int (skips all rows including header) to range (skips after header)
        chunk = pd.read_csv(input_file_path, skiprows=range(1, skip_rows), nrows=chunk_size, header=0)
        print(f"Retrying chunk {chunk_number}...")
        if process_chunk(chunk, chunk_number):
          print(f"Processed chunk {chunk_number} successfully on retry, sleeping for {sleep_time} seconds...")
          time.sleep(sleep_time)
        else:
          print(f"Failed to process chunk {chunk_number} again, logging for further review.")
          time.sleep(sleep_time)

# Example: Rerun incomplete chunks
rerun_incomplete_chunks(f"{output_dir}incomplete_chunks.json")

########################

# EDIT: Indicate which run # of incomplete/failed chunks
print(f"Incomplete chunks (3rd run): {incomplete_chunks}")
print(f"Failed chunks (3rd run): {failed_chunks}")

# Saving the list of incomplete/failed chunks
# EDIT: Indicate which run # of incomplete/failed chunks
with open(f"{output_dir}incomplete_chunks_run3.json", "w") as f:
    json.dump(incomplete_chunks, f)
with open(f"{output_dir}failed_chunks_run3.json", "w") as f:
    json.dump(failed_chunks, f)
