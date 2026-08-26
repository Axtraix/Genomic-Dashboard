import re
import requests
from typing import List, Dict, Set

# Standard HGNC gene symbol regex (2-8 uppercase letters/numbers starting with a letter)
GENE_SYMBOL_REGEX = r'\b[A-Z][A-Z0-9]{1,7}\b'

# Optional target reference set to eliminate common acronym false positives (e.g., DNA, RNA, USA, Phase)
KNOWN_GENES = {"PCSK9", "HTT", "HBB", "TTR", "ANGPTL3", "DNMT1", "CFTR", "BRCA1", "EGFR"}

def extract_gene_symbols(text: str, enforce_whitelist: bool = True) -> List[str]:
    """
    Extracts potential gene symbols from clinical trial titles and descriptions.
    """
    if not text:
        return []
    
    # Extract all uppercase candidate tokens
    candidates = set(re.findall(GENE_SYMBOL_REGEX, text))
    
    if enforce_whitelist:
        # Filter candidate tokens against known gene database
        found_genes = candidates.intersection(KNOWN_GENES)
    else:
        # Filter out common non-gene acronym noise
        noise_words = {"DNA", "RNA", "FDA", "NCT", "USA", "NIH", "OR", "AND", "NOT", "I", "II", "III", "IV"}
        found_genes = candidates - noise_words

    return sorted(list(found_genes))


def fetch_europe_pmc_patents(gene_symbol: str) -> Dict[str, any]:
    """
    Queries the Europe PMC REST API for patent literature matching 
    a given gene symbol and CRISPR/Gene Editing keywords.
    """
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    query_str = f'SRC:PAT AND ("CRISPR" OR "Gene Editing") AND "{gene_symbol}"'
    
    params = {
        "query": query_str,
        "format": "json",
        "pageSize": 50,
        "resultType": "core"  # Ensures full record details are returned
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        # Total number of matching patents found in Europe PMC
        total_hits = data.get("hitCount", 0)
        
        # List of actual patent records returned (up to 50)
        patent_list = data.get("resultList", {}).get("result", [])
        
        return {
            "total_count": total_hits,
            "retrieved_records_count": len(patent_list),
            "patents": patent_list
        }
    except requests.RequestException as e:
        print(f"API Error fetching patents for {gene_symbol}: {e}")
        return {"total_count": 0, "retrieved_records_count": 0, "patents": []}

# =====================================================================
# EXAMPLE INTEGRATION
# =====================================================================

sample_study_title = "A Phase 1 Study Evaluating Base Editing targeting PCSK9 and ANGPTL3 in Patients"
sample_description = "Investigating CRISPR-Cas9 therapeutic interventions for HTT gene silencing."

full_text = f"{sample_study_title} {sample_description}"

# 1. Parse Genes
detected_genes = extract_gene_symbols(full_text, enforce_whitelist=True)
print(f"Detected Gene Targets: {detected_genes}")

# 2. Fetch Patent Saturation per Gene
patent_results: Dict[str, int] = {}
for gene in detected_genes:
    patent_results[gene] = fetch_europe_pmc_patents(gene)

print("Patent Saturation Results:", patent_results)