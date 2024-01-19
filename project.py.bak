# necessary imports
from urllib.request import urlopen
import json
import collections
import os

# store the URL in url as parameter for urlopen
url = 'https://www.ebi.ac.uk/proteins/api/variation/P04637?format=json'

# store the response of URL
response = urlopen(url)

# storing the JSON response from url in data
data = json.loads(response.read())

# retriveing geneName, organismName and original sequence
geneName = data['geneName']
organismName = data['organismName']
sequence = data['sequence']

# initializing mutations dictionaries for storing them (known mutations is for the ones came from UniProt, generated mutations are the ones we did and finally allmutations is for combining them into 1 single dictionary)
allmutations = collections.defaultdict(dict)

# storing all aminoacids for comparisons
aminoacids = 'ACDEFGHIKLMNPQRSTVWY'

# iterating through the json for collecting every variant from UniProt
for features in data['features']:
    if features.get('begin') == features.get('end'):
        if (features.get('mutatedType') != None and (features.get('mutatedType') in aminoacids and features.get('wildType') in aminoacids)):
            variantSequence = sequence[:int(features.get('begin'))-1] + features.get('mutatedType') + sequence[int(features.get('begin')):]
            allmutations[variantSequence]['name'] = features.get('locations')[0].get('loc')
            allmutations[variantSequence]['description'] = features.get('descriptions')[0].get('value') if features.get('descriptions') != None else ""
            allmutations[variantSequence]['position'] = int(features.get('begin'))
            allmutations[variantSequence]['originalType'] = features.get('wildType')
            allmutations[variantSequence]['mutatedType'] = features.get('mutatedType')

# generating and storing mutations
count = 0
for original in range(0, len(sequence)):
    for mutated in range(0, len(aminoacids)):
        if sequence[original] == aminoacids[mutated]:
            continue
        else:
            newsequence = sequence[:original] + aminoacids[mutated] + sequence[original+1:]
            # for making sure not adding existing mutations again
            if newsequence not in allmutations.keys():
                allmutations[newsequence]['name'] = 'GeneratedMutation' + str(count)
                allmutations[newsequence]['description'] = ""
                allmutations[newsequence]['position'] = original+1
                allmutations[newsequence]['originalType'] = sequence[original]
                allmutations[newsequence]['mutatedType'] = aminoacids[mutated]
                count += 1

# generating files for doing queries on their respective algorithms
with open(os.path.join('prediction queries', 'tp53.fasta'), "w") as outfile:
    outfile.write('>'+geneName+' '+organismName)
    outfile.write('\n')
    outfile.write(sequence)
outfile.close()

with open(os.path.join('prediction queries', 'siftquery.txt'), "w") as outfile:
    for i in allmutations.keys():
        outfile.write(allmutations[i]['originalType']+str(allmutations[i]['position'])+allmutations[i]['mutatedType'])
        outfile.write('\n')
outfile.close()

with open(os.path.join('prediction queries', 'polyphenquery.txt'), "w") as outfile:
    for i in allmutations.keys():
        outfile.write('P04637'+' '+str(allmutations[i]['position'])+' '+allmutations[i]['originalType']+' '+allmutations[i]['mutatedType'])
        outfile.write('\n')
outfile.close()

# PolyPhen and SIFT matrix creations and storing them in their own respective dictionaries
siftpredictionmatrixtype = collections.defaultdict(dict)
siftpredictionmatrixvalue = collections.defaultdict(dict)
polyphen2predictionmatrixtype = collections.defaultdict(dict)
polyphen2predictionmatrixvalue = collections.defaultdict(dict)

with open(os.path.join('inputs', 'sift predictions.txt')) as f:
    for line in f:
        predictionvalues = line.split()
        for aminoacid in range(0, len(aminoacids)):
            if (float(predictionvalues[aminoacid+2])>0.05):
                predictionValType = 'tolerated'
            else:
                predictionValType = 'damaging'

            siftpredictionmatrixtype[predictionvalues[0]][aminoacids[aminoacid]] = predictionValType
            siftpredictionmatrixvalue[predictionvalues[0]][aminoacids[aminoacid]] = float(predictionvalues[aminoacid+2])

with open(os.path.join('inputs', 'polyphen predictions.txt')) as f:
   for line in f:
        predictionvalues = line.split()
        if (predictionvalues[9] == 'benign'):
            polyphen2predictionmatrixtype[predictionvalues[6]+predictionvalues[7]][predictionvalues[8]] = predictionvalues[9]
            polyphen2predictionmatrixvalue[predictionvalues[6]+predictionvalues[7]][predictionvalues[8]] = float(predictionvalues[10])
        else:
            polyphen2predictionmatrixtype[predictionvalues[6]+predictionvalues[7]][predictionvalues[8]] = predictionvalues[9] + ' ' + predictionvalues[10]
            polyphen2predictionmatrixvalue[predictionvalues[6]+predictionvalues[7]][predictionvalues[8]] = float(predictionvalues[11])

# adding prediction values for both algoritms into the dictionary
for i in allmutations.keys():
    allmutations[i]['siftpredictionType'] = siftpredictionmatrixtype[str(allmutations[i]['position'])+allmutations[i]['originalType']][allmutations[i]['mutatedType']]
    allmutations[i]['siftpredictionValue'] = siftpredictionmatrixvalue[str(allmutations[i]['position'])+allmutations[i]['originalType']][allmutations[i]['mutatedType']]
    allmutations[i]['polyphen2predictionType'] = polyphen2predictionmatrixtype[str(allmutations[i]['position'])+allmutations[i]['originalType']][allmutations[i]['mutatedType']]
    allmutations[i]['polyphen2predictionValue'] = polyphen2predictionmatrixvalue[str(allmutations[i]['position'])+allmutations[i]['originalType']][allmutations[i]['mutatedType']]

# creating json file for r
with open(os.path.join('outputs', 'dictionary.csv'), "w") as outfile:
    outfile.write("name,sequence,description,position,originalType,mutatedType,siftpredictionType,siftpredictionValue,polyphen2predictionType,polyphen2predictionValue")
outfile.close()

with open(os.path.join('outputs', 'dictionary.csv'), "a") as outfile:
    for i in allmutations.keys():
        outfile.write('\n')
        outfile.write(allmutations[i]['name']+','+i+','+allmutations[i]['description']+','+str(allmutations[i]['position'])+','+allmutations[i]['originalType']+','+allmutations[i]['mutatedType']+','+allmutations[i]['siftpredictionType']+','+str(allmutations[i]['siftpredictionValue'])+','+allmutations[i]['polyphen2predictionType']+','+str(allmutations[i]['polyphen2predictionValue']))
outfile.close()
