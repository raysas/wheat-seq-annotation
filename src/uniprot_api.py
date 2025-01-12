from Bio import ExPASy, SwissProt

def extract_chromosome_info(data):
    chromosome_info = []
    
    for entry in data:
        for element in entry:
            if 'Chromosome' in element:
                chromosome_info.append(element)
    
    return ''.join(chromosome_info)

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
        * Genomic Cross-references
    '''
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        # print(extract_chromosome_info(record.cross_references)) # success :D
        handle.close()
        
        annotations = {
            "Protein Name": record.description,
            "Genomic position": extract_chromosome_info(record.cross_references),
            "Organism": record.organism,
            "Entry Name": record.entry_name,
            "Gene Name": record.gene_name,
            "Keywords": record.keywords,
            "Comments": record.comments,
            "ID": uniprot_id
        }
        
        return annotations
    except Exception as e:
        return {"error": f"An error occurred: {e}"}
import sys

if __name__ == "__main__":
    #example: python uniprot_api.py A0A453GRS5
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
