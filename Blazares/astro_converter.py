
# Program by Paula Gálvez Molina on March -  2023.
# Defines functions for converting a single number or a whole file on time, flux
# etc.
# The path for the file must be an absolute (i.e. complete) path, not relative.
# The numpy and astropy libraries are required for the code to function.

# SUGGESTIONS: When converting time and flux, first convert and flux to magnitude
# and then convert the time column to keep more number of decimals.

import numpy as np
from astropy import units as u
from astropy.time import Time
import pandas as pd
from sklearn import preprocessing

def convert_MJD_to_MET(number_in_MJD_UTC):
    # This function converts the date provided in MJD [UTC] to Fermi's MET
    # by substracting the 0 time for the MET in MJD to the given date, to then
    # convert to and return the difference in seconds.

    # Define that 51910.00000000 MJD [UTC] = 0.000 MET
    MET_zero_in_MJD_UTC = Time(51910.00000000, format='mjd')

    # Assign the MJD [UTC] to the input time
    number_in_MJD_UTC = Time(number_in_MJD_UTC, format='mjd')

    # Take the difference, converting, formatting, and returing value
    difference =  number_in_MJD_UTC - MET_zero_in_MJD_UTC
    date_MET = f"{difference.to(u.second).value:0.03f}"
    return(date_MET)

def convert_File_MJD_to_MET(file_name):
    # Accepts a csv (or .dat) file where the first column is the time in MJD
    # [UTC] to be converted to Fermi's MET. The remaining columns of the file
    # remain unchanged, and returns a the same type of file it received with the
    # first column converted to MET under the original name +_converted.

    # Check the data file is either .csv or .dat
    if (file_name[-4:] != ".csv" and file_name[-4:] != ".dat"):
        print("File type is not supported, this function only takes .csv or .dat files.")
        return(0)

    # Load the input file and extract the time column
    data = np.genfromtxt(file_name, dtype=float)
    header = open(file_name).readline().rstrip()[1:]
    time_column = data[:, 0]

    # Convert the time column to MET using the vectorized function
    converted_column = np.vectorize(convert_MJD_to_MET)(time_column)

    # Combine the converted time column and the remaining columns of the input
    # file into a single NumPy array
    converted_data = np.concatenate((converted_column[:, np.newaxis], data[:, 1:]), axis=1)

    # Save the converted data to a new file
    np.savetxt(fname = file_name[:-4]+"_converted"+file_name[-4:],  X = converted_data, comments = '#', fmt="%s", delimiter="\t", newline="\n", header = header)

def convert_magnitude_to_flux(base_magnitude, base_flux, magnitude_to_convert):
    # Using Pogson's equation
    converted_flux = base_flux*100**((base_magnitude-magnitude_to_convert)/5)
    return(converted_flux)

def convert_File_magnitude_to_flux(base_magnitude, base_flux, file_name):
    # Accepts a csv (or .dat) file where the third column is the magnitude
    # to convert to flux. The remaining columns of the file
    # remain unchanged, and returns a the same type of file it received with the
    # third column converted to flux under the original name +_converted.

    # Check the data file is either .csv or .dat
    if (file_name[-4:] != ".csv" and file_name[-4:] != ".dat"):
        print("File type is not supported, this function only takes .csv or .dat files.")
        return(0)

    # Load the input file and extract the magnitude column
    data = np.genfromtxt(file_name, dtype=float)
    header = open(file_name).readline().rstrip()[1:]
    magnitude_column = data[:, 2]

    # Convert the flux column to magnitude using the vectorized function
    converted_column = np.vectorize(convert_magnitude_to_flux)(base_magnitude, base_flux, magnitude_column)

    # Combine the converted magnitude column and the remaining columns of the input
    # file into a single NumPy array
    converted_data = np.concatenate((data[:, :2], converted_column[:, np.newaxis], data[:, 3:]), axis=1)

    # Save the converted data to a new file
    # np.savetxt(fname = file_name[:-4]+"_converted"+file_name[-4:],  X = converted_data, comments = '#', fmt="%s", delimiter="\t", newline="\n", header = header)

def sort_File_by_column(file_name, column_number = 0, ascending = True,
separator = "\t"):
    # Accepts a csv (or .dat) with a commented header file to sort the entire
    # dataframe by the ascending or descending order of the indicated column.
    # Returns the same type of file it received with the
    # under the original name +_sorted.
    # By default the separator is a tab, but it might be changed to a comma or " ".
    # By default the given data frame will be sorted ascendingly by the values
    # of the first column.

    # Gets the header
    header = open(file_name).readline().rstrip()[1:]

    # Converts file to dataframe and sorts by specified column
    df = pd.read_csv(file_name, sep = separator, skiprows=1, header=None)
    first_col = df.iloc[:, column_number]
    df.sort_values(by=first_col.name, inplace = True)

    # Save the sorted data to a new file
    np.savetxt(fname = file_name[:-4]+"_sorted"+file_name[-4:],  X = df, comments = '#', fmt="%s", delimiter= "\t", newline="\n", header = header)

def Rodrigos_version_of_normalize_File_column(file_name, column_number = 2, separator = "\t"):
    # Accepts a csv (or .dat) with a commented header file to normalize a column
    # (by defect the third column) by dividing the value at that point by the
    # average of all.
    # Returns the same type of file it received with the
    # under the original name +_normalized.
    # By default the separator is a tab, but it might be changed to a comma or " ".
    # By default the given data frame will be sorted ascendingly by the values
    # of the first column.

    try:
        # Gets the header
        header = open(file_name).readline().rstrip()[1:]

        # Converts file to dataframe and sorts by specified column
        df = pd.read_csv(file_name, sep = separator, skiprows=1, header=None)

        # Check if the column_number is valid
        num_columns = df.shape[1]
        if column_number < 0 or column_number >= num_columns:
            raise ValueError(f"Invalid column_number. Must be in the range [0, {num_columns - 1}].")

        col = df.iloc[:, column_number]
        mean = np.mean(col)
        normalized_col = col/mean

        # Combine the normalized column and the remaining columns of the input
        # file into a single NumPy array
        converted_data = np.concatenate((df.iloc[:, :column_number], normalized_col[:, np.newaxis], df.iloc[:, column_number+1:]), axis=1)

        # Save the sorted data to a new file
        new_file_name = file_name[:-4] + "_normalized" + file_name[-4:]
        np.savetxt(fname = new_file_name,  X = converted_data, comments = '#', fmt="%s", delimiter= "\t", newline="\n", header = header)

        print(f"Normalization completed. The normalized file is saved as: {new_file_name}")

    # Error handling
    except FileNotFoundError:
            print(f"Error: File '{file_name}' not found.")
    except Exception as e:
            print(f"An error occurred: {e}")

def check_num_columns(file_name, separator='\t', skip_rows=1):
    try:
        # Read the file into a DataFrame
        df = pd.read_csv(file_name, sep=separator, skiprows=skip_rows)

        # Get the number of columns
        num_columns = df.shape[1]

        print(f"The file '{file_name}' has {num_columns} columns.")

    except FileNotFoundError:
        print(f"Error: File '{file_name}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

################################################################################
#                                   WORKSPACE                                  #
################################################################################
file_name = "/home/paula/Documents/Internship/AstroECFM/Blazares/PG1553+113/PG1553.dat"
# Rodrigos_version_of_normalize_File_column(file_name)
# sort_File_by_column(file_name)
# convert_File_magnitude_to_flux(1, 0.000000789, file_name)
# convert_File_MJD_to_MET(file_name)

################################################################################
#                                       LOG                                    #
################################################################################

# convert_File_MJD_to_MET(file_name)
# convert_File_magnitude_to_flux(1, 0.000000789, file_name)
#file_name = "/home/paula/Documents/Internship/Blazares/PG1553+113_multibands_normalized.dat"
#sort_File_by_column(file_name)
# Rodrigos_version_of_normalize_File_column(file_name)
#file1 = pd.read_csv("/home/paula/Documents/Internship/Blazares/PG1553+113_results_optical_converted_converted_sorted_normalized.dat", sep = "\t")
# file2 = pd.read_csv("/home/paula/Documents/Internship/Blazares/PG1553+113_results_45dias_gamma_normalized.dat", sep = "\t")
# combined_df = pd.concat([file1, file2[]], ignore_index=True)
# np.savetxt(fname = file_name[:-4]+"_normalized"+file_name[-4:],  X = converted_data, comments = '#', fmt="%s", delimiter= "\t", newline="\n", header = header)
