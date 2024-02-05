#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd


# In[3]:


m_path =  ['I:/LAB SAMPLES/1. LARGE-PD from LATAM/2. Brazil/Quantification Raw Data/01-2024/RAWDATA Ribeirao 8-22 stock tubes dsDNA 1_22_2024 10_50_44 AM.csv'] 
# input path(s) to the raw NanoDrop output .csv file(s) here
rawdata = []
# leave blank
for path in m_path:
    df1 = pd.read_csv(path,encoding='utf-16', delimiter='\t',usecols=['Nucleic Acid(ng/uL)'])
    df1['Consolidated Nano Raw'] = path
    rawdata.append(df1)
nano_rawdata = pd.concat(rawdata, ignore_index=True)
nano_rawlist = nano_rawdata['Nucleic Acid(ng/uL)'].tolist()
print(nano_rawlist)


# nanoen = list(enumerate(nano_rawlist))
# print(nanoen)

# In[5]:


source = input('What is the source of these readings? (type either SP or ST):')
#inputs should only equal 'SP' or 'ST'


# In[6]:


code_path = 'I:/LAB SAMPLES/1. LARGE-PD from LATAM/2. Brazil/Quantification Raw Data/01-2024/Ribeirao codes 8-22 1-22-24.csv'
# input path to .csv file for the list of codes exported by the NanoDrop
pcodes = pd.read_csv(code_path,usecols=[0,1])
pcodes_df = pd.DataFrame(pcodes)
pcodesen = list(enumerate(pcodes_df["codes"]))
pcodes_index = {code:volume for code, volume in zip(pcodes_df["codes"], pcodes_df["volumes"])}
poscode_index = {position:pdcode for position, pdcode in pcodesen}
print(poscode_index)
print(pcodes_index)


# In[9]:


nanoavg = []
# leave blank
index_map = {spot:round(nano,1) for spot, nano in nanoen}
skip_until = 0
cutoff = ((0,99,2),(99,299,5),(299,499,10),(499,999,15),(999,10000,20))
cutval = []
for spot, nano in nanoen:
    if spot < skip_until:
        continue
    if source == "SP" and index_map[spot] < 35:
        continue
    if source == "SP" and index_map[spot] > 95:
        continue
    for n1, n2, cut in cutoff:
        if (index_map[spot] > n1) and (index_map[spot] < n2 or index_map[spot] == n2):
            cutval = cut
            break
        else:
            continue
    if int(spot + 1) in index_map.keys():
        if abs(index_map[spot] - index_map[spot+1]) <= cutval:
            nanoavg.append((round((index_map[spot]+index_map[spot+1])/2,1)))
            skip_until = spot + 2
        else:
            combos = [(0,1),(0,2),(1,2),(0,3),(1,3),(2,3),
                      (0,4),(1,4),(2,4),(3,4),(0,5),(1,5),
                      (2,5),(3,5),(4,5),(0,6),(1,6),(2,6),
                      (3,6),(4,6),(5,6),(0,7),(1,7),(2,7),
                      (3,7),(4,7),(5,7),(6,7),(0,8),(1,8),
                      (2,8),(3,8),(4,8),(5,8),(6,8),(7,8)]
            for spot1,spot2 in combos:
                if source == "SP" and (index_map[spot + spot1] < 35 or index_map[spot + spot2] < 35 
                                       or index_map[spot + spot1] > 95 or index_map[spot +spot2] > 95):
                    continue
                if (spot + spot1) and (spot + spot2) in index_map.keys():
                    if (abs(index_map[spot + spot1]-index_map[spot + spot2])) <= cutval:
                        nanoavg.append((round((index_map[spot + spot1]+index_map[spot + spot2])/2,1)))
                        if spot1 > spot2:
                            skip_until = spot + spot1 + 1
                            break
                        elif spot2 > spot1:
                            skip_until = spot + spot2 + 1
                            break
print(nanoavg)


# In[43]:


index_map = {spot:round(nano,1) for spot, nano in nanoen}
print(index_map)


# In[59]:


def selectCutoff(concentration):
    if 0 < concentration <= 99:
        return 2
    elif 99 < concentration <= 299:
        return 5
    elif 299 < concentration <= 499:
        return 10
    elif 499 < concentration <= 999:
        return 15
    elif 999 < concentration:
        return 20
    else:
        input("Talk with the developer about concentration zero")


def checkIfWillDilute(source, concentration1, concentration2):
    if source == "SP" and (concentration1 < 35 or concentration2 < 35 or concentration1 > 95 or concentration2 > 95):
        return True
    return False


nanoavg = []
# leave blank

#Dictionary with the position in the colum and nano raw value with one decimal place
index_map = {spot:round(nano,1) for spot, nano in nanoen}
skip_until = 0
samplenum = 0
#cutval = []
for spot, nano in nanoen:
   # print(f"{source} {spot} {nano} {nanoavg}")
    
    #If fills any condition, ignore this line
    if spot < skip_until:
        #print("continue 1")
        continue
    if source == "SP" and (index_map[spot] < 35 or index_map[spot] > 95):
        #print("continue 2")
        continue
    
    #Select the cutoff based on the function
    cutval = selectCutoff(index_map[spot])
    
    while samplenum in poscode_index.keys():
        ID = poscode_index[samplenum]
        #print(f"{ID} - {pcodes_index[ID]}")
        if pcodes_index[ID] == 0:
            #print("Append zero")
            nanoavg.append(0)
            samplenum = samplenum + 1
        else:
            #print("break")
            break
            
            
    if int(spot + 1) in index_map.keys():
        ID = poscode_index[samplenum]
        if pcodes_index[ID] == 4:
            #print("4")
            nanoavg.append(index_map[spot])
            samplenum = samplenum + 1
            skip_until = spot + 1
        elif pcodes_index[ID] < 10:
            #print("<10")
            nanoavg.append((round((index_map[spot]+index_map[spot+1])/2,1)))
            samplenum = samplenum + 1
            skip_until = spot + 2
        elif pcodes_index[ID] == 10:
            #print("=10")
            combos = [(0,2),(1,2)]
            for spot1,spot2 in combos:
                if (spot + spot1) and (spot + spot2) in index_map.keys():
                    if (abs(index_map[spot + spot1]-index_map[spot + spot2])) <= 100:
                        nanoavg.append((round((index_map[spot + spot1]+index_map[spot + spot2])/2,1)))
                        skip_until = spot + spot2 + 1
                        samplenum = samplenum + 1
                        break
        elif abs(index_map[spot] - index_map[spot+1]) <= cutval:
            #print("Cutval")
            nanoavg.append((round((index_map[spot]+index_map[spot+1])/2,1)))
            samplenum = samplenum + 1
            skip_until = spot + 2
        else:
            #minimum indexes
            spot1=0
            spot2=1

            #Set dictionary with already compared values
            dictLooked={}
            while spot2 in index_map:
                #Set the first key
                if spot1 not in dictLooked:
                    dictLooked[spot1]= []
                
                #Check if I alread compared spot1 and spot2
                if spot2 not in dictLooked[spot1]:
                    dictLooked[spot1].append(spot2)
                else:
                    if checkIfWillDilute(source, index_map[spot + spot1], index_map[spot + spot2]):
                        #print("Dilution")
                        continue
                    if (spot + spot1) and (spot + spot2) in index_map.keys():
                        #print(f"{spot} {spot1} {spot2}")
                        
                        cutvalFirst = selectCutoff(index_map[spot + spot1])
                        cutvalSecond = selectCutoff(index_map[spot + spot2])
                        
                        cutval = cutvalSecond
                        if cutvalFirst > cutvalSecond :
                            cutval=cutvalFirst
                        
                        #print(f"abs({index_map[spot + spot1]}-{index_map[spot + spot2]}) <= {cutval}")
                        if (abs(index_map[spot + spot1]-index_map[spot + spot2])) <= cutval:
                            #print(f"Lower than cutoff {cutval} with {spot} {spot1} {spot2}")
                            nanoavg.append((round((index_map[spot + spot1]+index_map[spot + spot2])/2,1)))
                            skip_until = spot + spot2 + 1
                            samplenum = samplenum + 1
                            break
                    
                    
                    if spot2-spot1 >= 2:
                        spot1 = spot1+1
                    else:
                        spot2=spot2+1
                        spot1=0
    #print(nanoavg)
    #input()
            
while samplenum in poscode_index.keys():
        ID = poscode_index[samplenum]
        if pcodes_index[ID] == 0:
            nanoavg.append(0)
            samplenum = samplenum + 1
            continue
        else:
            break
print(nanoavg)


# In[62]:


code_path = 'I:/LAB SAMPLES/1. LARGE-PD from LATAM/2. Brazil/Quantification Raw Data/01-2024/Ribeirao codes 8-22 1-22-24.csv'
# input path to .csv file for the list of codes exported by the NanoDrop
pcodes= pd.read_csv(code_path,usecols=['codes','volumes'])

codes1 = np.array([pcodes])
pcodes['Nanodrop Average [ng/uL]'] = nanoavg
print(pcodes)


# In[63]:


df2 = pd.DataFrame(pcodes)
df2.to_csv('I:/LAB SAMPLES/1. LARGE-PD from LATAM/2. Brazil/Quantification Raw Data/01-2024/Ribeirao nano avg codes 8-22 1-22-24.csv')
# input the path to where the output file should go including the name of the output file

