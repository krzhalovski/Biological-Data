import requests
import time
import sys
import pandas as pd


# Example: python3 1.1_pfam_domains.py https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/PF00483/ reviewed_domain_positions.csv

URL = sys.argv[1]
OUT = sys.argv[2]

attempts = 0

results = []
i = 0

while True:
    if attempts > 50:
        break
    try:
        response = requests.get(URL)
    except:
        print(f"Connection problem, waiting 5 minutes {URL}")
        time.sleep(5*60)
        attempts += 1
        continue
    
    if response.status_code != 200:
        print(f"Response error: {response.status_code}, waiting 5 minutes {URL}")
        time.sleep(5*60)
        attempts += 1
        continue

    info = response.json()

    URL = info['next']

    for result in info['results']:
        accession = result['metadata']['accession']

        for entry in result['entries']:
            for location in entry['entry_protein_locations']:
                if location['model'] != 'PF00483':
                    continue

                for fragment in location['fragments']:
                    start = fragment['start']
                    end = fragment['end']

                    results.append({'accession': accession, 'start': start, 'end': end})

    if i%10==0:
        print(f"Pages parsed: {i}")

    i += 1
    time.sleep(1)

    if URL == None:
        break

print(f"Total pages parsed: {i}")
df = pd.DataFrame(results)
df.to_csv(OUT)