import requests
from functools import lru_cache
from flask import Flask, jsonify, send_file
from flask_cors import CORS
from typing import List, Dict, Any

app = Flask(__name__, static_folder="static")
CORS(app)  # Enables cross-origin requests for local development

KEYWORD_TO_GENE = {
    "INSULIN": "INS",
    "GLP-1": "GLP1R",
    "GLP1": "GLP1R",
    "GLUCAGON": "GCGR",
    "PCSK9": "PCSK9",
    "CFTR": "CFTR",
    "EGFR": "EGFR",
    "BRCA1": "BRCA1",
    "TTR": "TTR"
}

@lru_cache(maxsize=256)
def fetch_patent_count(gene_symbol: str) -> int:
    if gene_symbol == "None":
        return 0
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    query_str = f'SRC:PAT AND ("CRISPR" OR "Gene Editing" OR "Therapy") AND "{gene_symbol}"'
    params = {"query": query_str, "format": "json", "pageSize": 1}
    try:
        response = requests.get(url, params=params, timeout=5)
        response.raise_for_status()
        return response.json().get("hitCount", 0)
    except Exception:
        return 0

def extract_gene(text: str) -> str:
    if not text:
        return "None"
    text_upper = text.upper()
    for kw, gene in KEYWORD_TO_GENE.items():
        if kw in text_upper:
            return gene
    return "None"

def fetch_clinical_trials_data(limit: int = 200) -> List[Dict[str, Any]]:
    # Fetching up to 200 items in a single reliable request to avoid pagination rate-limits
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.cond": "Diabetes",
        "filter.overallStatus": "RECRUITING",
        "pageSize": min(limit, 100)
    }

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
    except Exception as e:
        print(f"[API ERROR] ClinicalTrials fetch failed: {e}")
        return []

    studies = data.get('studies', [])
    records = []

    for s in studies:
        protocol = s.get('protocolSection', {})
        id_mod = protocol.get('identificationModule', {})
        status_mod = protocol.get('statusModule', {})
        desc_mod = protocol.get('descriptionModule', {})
        cond_mod = protocol.get('conditionsModule', {})
        sponsor_mod = protocol.get('sponsorCollaboratorsModule', {})

        title = id_mod.get('briefTitle', 'N/A')
        summary = desc_mod.get('briefSummary', '')
        detected_gene = extract_gene(f"{title} {summary}")

        records.append({
            "nct_id": id_mod.get('nctId', 'N/A'),
            "title": title,
            "status": status_mod.get('overallStatus', 'N/A'),
            "primary_gene": detected_gene,
            "conditions": ", ".join(cond_mod.get('conditions', [])[:2]),
            "sponsor": sponsor_mod.get('leadSponsor', {}).get('name', 'N/A')
        })

    if not records:
        return []

    # Calculate metrics
    gene_counts = {}
    for r in records:
        g = r['primary_gene']
        gene_counts[g] = gene_counts.get(g, 0) + 1

    patent_counts = {g: fetch_patent_count(g) for g in gene_counts}

    for record in records:
        target = record['primary_gene']
        patents = patent_counts.get(target, 0)
        trials = gene_counts.get(target, 1)
        idx = round(patents / trials, 2) if target != "None" else 0.0

        if target == "None":
            status_label = "Unassigned"
        elif idx > 50:
            status_label = f"Saturated ({idx})"
        else:
            status_label = f"Unmet Need ({idx})"

        record["patent_count"] = patents
        record["saturation_index"] = idx
        record["ip_index_label"] = status_label

    return records

@app.route("/")
def index():
    return send_file("crispr_dashboard.html")

@app.route("/api/studies")
def api_studies():
    data = fetch_clinical_trials_data(limit=100)
    return jsonify(data)

if __name__ == "__main__":
    print("Starting Analytics Flask Server on http://127.0.0.1:5000 ...")
    app.run(debug=True, port=5000)