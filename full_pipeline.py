import re
import requests
import pandas as pd
from flask import Flask, jsonify, send_file
from typing import List

app = Flask(__name__, static_folder="static")

GENE_SYMBOL_REGEX = r'\b[A-Z][A-Z0-9]{1,7}\b'
KNOWN_GENES = {"PCSK9", "HTT", "HBB", "TTR", "ANGPTL3", "DNMT1", "CFTR", "BRCA1", "EGFR", "INS", "INSR"}


def extract_gene_symbols(text: str) -> List[str]:
    if not text:
        return []
    candidates = set(re.findall(GENE_SYMBOL_REGEX, text))
    found_genes = candidates.intersection(KNOWN_GENES)
    return sorted(list(found_genes))


def fetch_patent_count(gene_symbol: str) -> int:
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    query_str = f'SRC:PAT AND ("CRISPR" OR "Gene Editing") AND "{gene_symbol}"'
    params = {"query": query_str, "format": "json", "pageSize": 1}
    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data.get("hitCount", 0)
    except requests.RequestException:
        return 0


def fetch_clinical_trials_data(limit: int = 50) -> List[dict]:
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.cond": "Diabetes",
        "filter.overallStatus": "RECRUITING",
        "pageSize": limit
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
    except requests.RequestException:
        return []

    raw_studies = data.get('studies', [])
    records = []

    for s in raw_studies:
        protocol = s.get('protocolSection', {})
        id_mod = protocol.get('identificationModule', {})
        status_mod = protocol.get('statusModule', {})
        desc_mod = protocol.get('descriptionModule', {})
        cond_mod = protocol.get('conditionsModule', {})
        sponsor_mod = protocol.get('sponsorCollaboratorsModule', {})

        title = id_mod.get('briefTitle', 'N/A')
        summary = desc_mod.get('briefSummary', '')
        genes = extract_gene_symbols(f"{title} {summary}")

        records.append({
            "nct_id": id_mod.get('nctId', 'N/A'),
            "title": title,
            "status": status_mod.get('overallStatus', 'N/A'),
            "genes": genes,
            "conditions": ", ".join(cond_mod.get('conditions', [])),
            "sponsor": sponsor_mod.get('leadSponsor', {}).get('name', 'N/A')
        })

    df = pd.DataFrame(records)
    if df.empty:
        return []

    # Calculate IP Saturation Index per target gene
    gene_counts = {}
    for gene_list in df['genes']:
        for gene in gene_list:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1

    patent_cache = {gene: fetch_patent_count(gene) for gene in gene_counts}
    saturation_cache = {
        gene: round(patent_cache[gene] / gene_counts[gene], 2)
        for gene in gene_counts
    }

    processed_records = []
    for record in records:
        genes = record['genes']
        if genes:
            target = genes[0]
            patents = patent_cache.get(target, 0)
            idx = saturation_cache.get(target, 0.0)
            status_label = f"Saturated (Index: {idx})" if idx > 50 else f"High Unmet Need (Index: {idx})"
        else:
            target, patents, status_label = "None", "N/A", "N/A"

        record["primary_gene"] = target
        record["patent_count"] = patents
        record["ip_index_label"] = status_label
        processed_records.append(record)

    return processed_records


@app.route("/")
def index():
    return send_file("index.html")


@app.route("/api/studies")
def api_studies():
    data = fetch_clinical_trials_data(limit=50)
    return jsonify(data)


if __name__ == "__main__":
    app.run(debug=True)