from Bio import ExPASy, SwissProt

def fetch_uniprot_annotation(uniprot_id:str)->dict:
    '''
    params:
        uniprot_id: str  
    returns:    
        annotations: dict of keys:  
        * Entry Name  
        * Protein Name  
        * Organism  
        * Gene Name  
        * Keywords  
        * Comments  
    '''
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        handle.close()
        
        annotations = {
            "ID": uniprot_id,
            "Entry Name": record.entry_name,
            "Protein Name": record.description,
            "Organism": record.organism,
            "Gene Name": record.gene_name,
            "Keywords": record.keywords,
            "Comments": record.comments,
        }
        
        return annotations
    except Exception as e:
        return {"error": f"An error occurred: {e}"}

import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python uniprot_api.py <UniProt ID>")
        sys.exit(1)

    uniprot_id = sys.argv[1].strip()
    annotation = fetch_uniprot_annotation(uniprot_id)

    if "error" in annotation:
        print(annotation["error"])
    else:
        for key, value in annotation.items():
            print(f"{key}: {value}")
