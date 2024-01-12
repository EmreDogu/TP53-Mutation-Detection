from urllib.request import urlopen
import json
import collections

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
            # for making sure not adding same mutations again
            if newsequence not in allmutations.keys():
                allmutations[newsequence]['name'] = 'GeneratedMutation' + str(count)
                allmutations[newsequence]['description'] = ""
                allmutations[newsequence]['position'] = original+1
                allmutations[newsequence]['originalType'] = sequence[original]
                allmutations[newsequence]['mutatedType'] = aminoacids[mutated]
                count += 1

#PolyPhen and SIFT matrix creations and adding prediction values to combined mutations dictionary
"""idk = []

            if features.get('predictions') != None:
                for predictions in features.get('predictions'):
                    if (predictions.get('predAlgorithmNameType') == 'SIFT' and predictions.get('version') != None) or predictions.get('predAlgorithmNameType') == 'PolyPhen':
                        idk.append([predictions.get('predAlgorithmNameType')  + ' ' + predictions.get('version'), predictions.get('predictionValType'), predictions.get('score')])
                        knownmutations[variantSequence]['predictions'] = idk
                    else:
                        idk.append([predictions.get('predAlgorithmNameType'), predictions.get('predictionValType'), predictions.get('score')])
                        knownmutations[variantSequence]['predictions'] = idk"""

#create csv files for r
