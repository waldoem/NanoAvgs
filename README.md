The JupyterLab file 'NanoAvgs.ipynb' was designed for the purpose of converting csv files from the NanoDrop One to a list of sample IDs and averages.
To operate use the following procedure:
  1. Save the .csv file from the NanoDrop One to the Mata Lab server in the appropriate site's 'Quantification Raw Data' folder
  2. Create a .csv file with a list of codes used with the header 'codes' and copy ONLY codes that were diluted for SPs (and therefore measured on the NanoDrop) and MAKE SURE the codes are in order of measurement. It will be helpful to copy the list directly from the dilution prep file for SPs
  4. Copy the path to the .csv file to the definition of m_path in the second cell
  5. Change all "\" to "/" in the path
  6. Copy the path for the .csv file with the codes in it to the 6th cell in the definition for 'code_path'
  7. Change all "\" to "/" in the path
  8. Change the location and name of the file for the output (this is VERY important so the previous file doesn't get overwritten)
  9. Run the code! make sure that the number of averages caluculated matches the number of codes in the code list (There should be an error when combining is attempted if they are not the same index length)
  10. Use a vlookup to copy Nano avgs to Master List and dilution files (when applicable)

Note: This code will not work if you take readings out of order of the sample list you import. 
