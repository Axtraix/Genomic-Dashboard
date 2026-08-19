import os
import re
import math
import requests
import pandas as pd
import numpy as np
import plotly.express as px
from typing import List, Optional, Dict
from pydantic import BaseModel, Field



class TrialRecord(BaseModel):
    nct_id: str = Field(..., description="Unique trial ID")
    title: str
    official_title: Optional[str] = "N/A"
    status: str
    phase: str
    start_date: Optional[str] = "N/A"
    sponsor: str
    sponsor_type: str
    conditions: List[str] = []
    interventions: List[str] = []

# Target gene database with prevalence scores (1-10 scale)
TARGET_GENES = {
    "PCSK9": {"disease": "Hypercholesterolemia", "prevalence": 9.5},
    "HTT":   {"disease": "Huntington's Disease", "prevalence": 4.0},
    "HBB":   {"disease": "Sickle Cell / Thalassemia", "prevalence": 8.0},
    "TTR":   {"disease": "ATTR Amyloidosis", "prevalence": 5.5},
    "ANGPTL3": {"disease": "Dyslipidemia", "prevalence": 8.5},
    "DNMT1": {"disease": "Hereditary Sensory Neuropathy", "prevalence": 2.5},
    "CFTR":  {"disease": "Cystic Fibrosis", "prevalence": 7.0}
}




def parse_clinical_study(study_data: dict) -> TrialRecord:
    """Parses raw JSON study objects from ClinicalTrials.gov v2 API."""
    protocol = study_data.get("protocolSection", {})
    
    id_info = protocol.get("identificationModule", {})
    status_info = protocol.get("statusModule", {})
    sponsor_info = protocol.get("sponsorCollaboratorsModule", {})
    design_info = protocol.get("designModule", {})
    conditions_info = protocol.get("conditionsModule", {})
    interventions_info = protocol.get("armsInterventionsModule", {})
    
    phase_list = design_info.get("phases", ["Not Specified"])
    phase_str = ", ".join(phase_list) if isinstance(phase_list, list) else str(phase_list)
    
    raw_interventions = interventions_info.get("interventions", [])
    intervention_names = [item.get("name", "") for item in raw_interventions if item.get("name")]
    
    return TrialRecord(
        nct_id=id_info.get("nctId", "N/A"),
        title=id_info.get("briefTitle", "N/A"),
        official_title=id_info.get("officialTitle", "N/A"),
        status=status_info.get("overallStatus", "N/A"),
        phase=phase_str,
        start_date=status_info.get("startDateStruct", {}).get("date", "N/A"),
        sponsor=sponsor_info.get("leadSponsor", {}).get("name", "N/A"),
        sponsor_type=sponsor_info.get("leadSponsor", {}).get("class", "N/A"),
        conditions=conditions_info.get("conditions", []),
        interventions=intervention_names
    )


def fetch_trials(search_keywords: List[str], max_pages: int = 2) -> pd.DataFrame:
    """Fetches gene-editing studies from ClinicalTrials.gov."""
    print("[1/4] Pulling trial data from ClinicalTrials.gov API...")
    url = "https://clinicaltrials.gov/api/v2/studies"
    results = []
    
    for query in search_keywords:
        token = None
        for _ in range(max_pages):
            params = {"query.term": query, "pageSize": 50, "format": "json"}
            if token:
                params["pageToken"] = token
                
            res = requests.get(url, params=params)
            if res.status_code != 200:
                print(f"Warning: Failed request for term '{query}'")
                break
                
            payload = res.json()
            for study in payload.get("studies", []):
                try:
                    record = parse_clinical_study(study)
                    results.append(record.model_dump())
                except Exception:
                    continue
                    
            token = payload.get("nextPageToken")
            if not token:
                break
                
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.drop_duplicates(subset=["nct_id"]).reset_index(drop=True)
    return df


def fetch_patents(genes: List[str]) -> Dict[str, int]:
    """Queries Europe PMC REST API for patent hit counts per gene target."""
    print("[2/4] Querying patent counts via Europe PMC API...")
    counts = {}
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    
    for gene in genes:
        query_str = f'SRC:PAT AND ("CRISPR" OR "Gene Editing") AND "{gene}"'
        params = {"query": query_str, "format": "json", "pageSize": 1}
        
        try:
            res = requests.get(base_url, params=params, timeout=10)
            if res.status_code == 200:
                data = res.json()
                counts[gene] = data.get("hitCount", 0)
            else:
                counts[gene] = 0
        except Exception:
            counts[gene] = 0
            
    return counts




def tag_genes_in_trials(df: pd.DataFrame, gene_list: List[str]) -> pd.DataFrame:
    """Extracts target genes from clinical trial titles and descriptions."""
    def match_gene(row):
        combined_text = f"{row['title']} {row['official_title']} {' '.join(row['conditions'])} {' '.join(row['interventions'])}"
        found = []
        for g in gene_list:
            if re.search(r'\b' + re.escape(g) + r'\b', combined_text, re.IGNORECASE):
                found.append(g)
        return found if found else ["Uncategorized"]

    df["detected_genes"] = df.apply(match_gene, axis=1)
    return df


def build_opportunity_matrix(df: pd.DataFrame, patent_counts: Dict[str, int]) -> pd.DataFrame:
    """Calculates the Opportunity Index score for each gene target."""
    print("[3/4] Building target gene analytics matrix...")
    
    exploded = df.explode("detected_genes")
    gene_trial_counts = exploded["detected_genes"].value_counts().to_dict()
    
    rows = []
    for gene, info in TARGET_GENES.items():
        trials_count = gene_trial_counts.get(gene, 0)
        patents_count = patent_counts.get(gene, 0)
        prev_score = info["prevalence"]
        
        # Opportunity Score = Prevalence / ((Trials + 1) * log10(Patents + 10))
        denominator = (trials_count + 1) * math.log10(patents_count + 10)
        score = round(prev_score / denominator, 2)
        
        rows.append({
            "Gene Target": gene,
            "Indication": info["disease"],
            "Active Trials": trials_count,
            "Patent Count": patents_count,
            "Prevalence Score": prev_score,
            "Opportunity Index": score
        })
        
    return pd.DataFrame(rows)




def render_dashboard(trials_df: pd.DataFrame, analytics_df: pd.DataFrame):
    """Renders interactive Plotly charts and outputs 'crispr_dashboard.html'."""
    print("[4/4] Rendering charts and generating HTML dashboard...")
    
    # Chart 1: Trials by Phase
    fig_phase = px.bar(
        trials_df["phase"].value_counts().reset_index(),
        x="phase",
        y="count",
        labels={"phase": "Trial Phase", "count": "Total Studies"},
        title="Active Clinical Trials by Phase",
        template="plotly_dark"
    )
    
    # Chart 2: Patent Saturation vs Clinical Activity
    fig_scatter = px.scatter(
        analytics_df,
        x="Patent Count",
        y="Active Trials",
        size="Prevalence Score",
        color="Opportunity Index",
        hover_name="Gene Target",
        text="Gene Target",
        title="IP Density vs Clinical Activity",
        template="plotly_dark"
    )
    fig_scatter.update_traces(textposition="top center")

    phase_html = fig_phase.to_html(full_html=False, include_plotlyjs="cdn")
    scatter_html = fig_scatter.to_html(full_html=False, include_plotlyjs="cdn")
    
    table_analytics = analytics_df.to_html(index=False, classes="table table-dark table-striped")
    table_trials = trials_df[["nct_id", "title", "phase", "status", "sponsor"]].head(10).to_html(index=False, classes="table table-dark table-sm")

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>CRISPR / Gene Therapy IP Dashboard</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{ background-color: #0d1117; color: #c9d1d9; font-family: sans-serif; padding: 20px; }}
            .card {{ background-color: #161b22; border: 1px solid #30363d; margin-bottom: 20px; padding: 15px; border-radius: 8px; }}
            .stat-number {{ font-size: 2rem; font-weight: bold; color: #58a6ff; }}
            .stat-label {{ font-size: 0.85rem; color: #8b949e; text-transform: uppercase; }}
        </style>
    </head>
    <body>
        <div class="container-fluid">
            <h2 class="mb-1">CRISPR & Gene Therapy Intelligence Dashboard</h2>
            <p class="text-secondary">Clinical trial tracking & target gene IP saturation analysis</p>
            <hr class="border-secondary">

            <!-- KPI Row -->
            <div class="row">
                <div class="col-md-3">
                    <div class="card text-center">
                        <div class="stat-number">{len(trials_df)}</div>
                        <div class="stat-label">Total Studies Analyzed</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="card text-center">
                        <div class="stat-number">{len(TARGET_GENES)}</div>
                        <div class="stat-label">Target Genes Tracked</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="card text-center">
                        <div class="stat-number">{analytics_df['Patent Count'].sum()}</div>
                        <div class="stat-label">Indexed Patents</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="card text-center">
                        <div class="stat-number">{analytics_df.loc[analytics_df['Opportunity Index'].idxmax()]['Gene Target']}</div>
                        <div class="stat-label">Top Market Opportunity</div>
                    </div>
                </div>
            </div>

            <!-- Opportunity Matrix -->
            <div class="card">
                <h5>Target Gene Opportunity Matrix</h5>
                <div class="table-responsive">{table_analytics}</div>
            </div>

            <!-- Plotly Charts -->
            <div class="row">
                <div class="col-md-6"><div class="card">{phase_html}</div></div>
                <div class="col-md-6"><div class="card">{scatter_html}</div></div>
            </div>

            <!-- Recent Ingested Trials -->
            <div class="card">
                <h5>Recent Ingested Studies</h5>
                <div class="table-responsive">{table_trials}</div>
            </div>
        </div>
    </body>
    </html>
    """

    with open("crispr_dashboard.html", "w", encoding="utf-8") as f:
        f.write(html_content)
        
    print("\n[Success] Created 'crispr_dashboard.html'. Double-click the file to open it in your web browser.")



if __name__ == "__main__":
    search_terms = ["CRISPR", "Gene Editing", "Base Editing"]
    genes = list(TARGET_GENES.keys())
    
    # Step 1: Ingest Data
    df_trials = fetch_trials(search_terms, max_pages=2)
    patent_counts = fetch_patents(genes)
    
    # Step 2: Process & Extract
    df_trials = tag_genes_in_trials(df_trials, genes)
    df_analytics = build_opportunity_matrix(df_trials, patent_counts)
    
    # Step 3: Export CSVs
    df_trials.to_csv("trials_dataset.csv", index=False)
    df_analytics.to_csv("opportunity_matrix.csv", index=False)
    
    # Step 4: Generate HTML Dashboard
    render_dashboard(df_trials, df_analytics)