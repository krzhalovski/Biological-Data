import gzip
import json
from lxml import etree

annotations = {}

with gzip.open('../data/SwissProt/uniprot_sprot.xml.gz') as sp:
    single_entry = ""
    
    for line in sp:
        decoded = line.decode()
        if decoded.startswith('<entry'):
            single_entry = ""
        single_entry += decoded
        if decoded.startswith('</entry>'):
            root = etree.fromstring(single_entry)
            namespaces = root.nsmap
            entry_annotations = []

            for dbReference in root.findall('dbReference', namespaces):
                if dbReference.get('type') == 'GO':
                    entry_annotations.append(dbReference.get('id'))

            if entry_annotations:
                annotations[root.find("accession", namespaces).text] = entry_annotations
                single_entry = ""

with open("../data/SwissProt/annotations.json", "w") as f:
    json.dump(annotations, f)

breakpoint()


