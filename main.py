import csv
import numpy as np
import pandas as pd
import statistics
import os
os.getcwd() # provides the full file path for the current working directory
os.listdir() # provides a list of all files in the current working directory
global final_name

folders = "/Users/margauxvasilescu"
#This function matches the name of the weights with the pdb
def get_file_names_with_strings(str_list):

    full_list = os.listdir(f"{folders}/ConvertTensorWeights/done_matched_pdb")
    for word in full_list:
        #print("this is word", word)
        if not word.startswith('.') and str_list in word:
            final_name = word
            #print("this is final_name", final_name)

            return final_name

#this function fixes the texts.
def text_num_split(item):
    for index, letter in enumerate(item, 0):
        if letter.isdigit():
            return [item[:index],item[index:]]

pos_orig = {'Q15306': 428, 'Q8TDN2': 537, 'Q9P2S6': 618, 'Q6P2D8': 678, 'P51168': 400, 'Q6VMQ6': 214, 'O75533': 902, 'P14679': 79, 'P04637': 205, 'Q8N127': 318}

#this with statements creates a csv file to put the data that matches the pdbs that have highest weights that are close the pos of origin.
with open(f'{folders}/FilterWeights/csv_file.csv', 'w') as f:
    # create the csv writer
    writer = csv.writer(f)

    header = ['row', 'pos_orig', 'alpha_num','top_weight', 'weight_value', 'filename']
    writer.writerow(header)


    n = 0
    directory = f'{folders}/ConvertTensorWeights/new_weights_folder3'
    check_name = "xccccccccccc"
    list_checked_name = ["xccccccccccc"]
    seven_val_list = []
    eight_val_list = []
    nine_val_list = []
    alpha_list = []


    #this for loop goes through the pdb files with the weights to find the files with the highest weights that are closest to the pos of orig.
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)

        #this if statement makes sure that the file is a real pdb and doesn't starts with a period.
        if not filename.startswith('.') and os.path.isfile(f):
            pdb_name = os.path.basename(f).split("_")[3]
            protein_file = os.path.join(f'{folders}/ConvertTensorWeights/new_weights_folder3/{filename}')
            #print("this is protein_file", protein_file)
            #This finds the matching pos_orig in the file.
            pos_orig_no = pos_orig[pdb_name]
            pos_orig_no = int(pos_orig_no)

            outfile = open(protein_file, "r")
            data = outfile.readlines()
            new_data = []
            new_data3 = []
            alpha_list = []

            total_weight_data = []
            n = n + 1

            for line in data:
                if "ATOM" in line:
                   #print("line with atom", line)
                   line_new = line.split()
                   #print("new atom protein line number", line_new[5])
                   #print("new atom protein line", line)
                   new_data.append(line)
            new_data2 = np.array(new_data)

            new_data_weights = []


            for i in new_data2:
                #print("this is line x", i)
                line_new = i.split()
                if len(line_new) != 12:
                        if len(line_new[4]) != 1:
                            #print("spot 5", line_new[4])
                            insert_split = text_num_split(line_new[4])
                            line_new.insert(4, insert_split[0])
                            line_new[5] = insert_split[1]

                        spot_list = 0
                        for string in line_new:
                            t = string.split('.')
                            #print("this is t", t)
                            if len(t)>2:
                                splt_char = "-"
                                K = 2
                                # Printing original string
                                #print("The original string is : " + str(string))

                                # Split string on Kth Occurrence of Character
                                # using split() + join()
                                temp = string.split(splt_char)
                                res = splt_char.join(temp[:K]), splt_char.join(temp[K:])

                                string = res
                                line_new.insert(spot_list, res[0])
                                second_spot = spot_list + 1
                                line_new[second_spot] = res[1]


                            spot_list = spot_list + 1


                new_data3.append(line_new)
                six_value = int(line_new[5])
                weight_value = line_new[10]
                weight_value = float(weight_value)
                num_row = six_value
                num_row = int(num_row)
                total_weight_data.append(weight_value)





            top_weight = max(total_weight_data)
            top_weight = 0.90 * top_weight
            new_data4 = np.array(new_data3)

            for i in new_data4:
                seven_val = i[6]
                eight_val = i[7]
                nine_val =  i[8]



                if str(i[5]) == str(pos_orig_no):
                    seven_val_list.append(float(seven_val))
                    eight_val_list.append(float(eight_val))
                    nine_val_list.append(float(nine_val))
                    seven_ave = sum(seven_val_list)/len(seven_val_list)
                    eight_ave = sum(eight_val_list)/len(eight_val_list)
                    nine_ave = sum(nine_val_list)/len(nine_val_list)



            times = 0
            for i in new_data4:
                    times = times + 1
                    weight_value = i[10]
                    weight_value = float(weight_value)
                    #print("this is weight_value and top weight", weight_value, top_weight)
                    #print("this is six_value and v", six_value)

                    num_row = six_value
                    num_row = int(num_row)
                    total_weight_data.append(weight_value)
                    seven_no = float(i[6])
                    eight_no = float(i[7])
                    nine_no =  float(i[8])
                    alpha_num = ((seven_no - seven_ave)**2 + (eight_no - eight_ave)**2 + (nine_no - nine_ave)**2)**0.5
                    alpha_list.append(alpha_num)


            max_alpha = max(alpha_list)
            min_alpha = min(alpha_list)

            range_alpha = max_alpha - min_alpha
            max_range = (0.10 * range_alpha) + min(alpha_list)


            times = 0
            for i in new_data4:
                    times = times + 1
                    seven_no = float(i[6])
                    eight_no = float(i[7])
                    nine_no =  float(i[8])
                    alpha_num = ((seven_no - seven_ave)**2 + (eight_no - eight_ave)**2 + (nine_no - nine_ave)**2)**0.5
                    weight_value = i[10]
                    weight_value = float(weight_value)

                    if  (alpha_num <= max_range) and weight_value >= top_weight:

                        if pdb_name not in list_checked_name:
                            print("this is i in new_data4", i)
                            print("match",  pdb_name, filename)
                            print("list checked name", list_checked_name)
                            row = [i[5], pos_orig_no, alpha_num, top_weight, weight_value, filename]
                            writer.writerow(row)
                            print("alpha num", alpha_num)
                            print("max alpha", max_alpha)
                            print("min alpha", min_alpha)
                            print("max range", max_range)
                            print("weight value", weight_value)
                            print("top weight", top_weight)

                        if pdb_name not in list_checked_name:

                            list_checked_name.append(pdb_name)



