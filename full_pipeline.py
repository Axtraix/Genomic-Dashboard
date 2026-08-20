import os
import re
import math
import logging
import requests
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from typing import List, Optional, Dict
from pydantic import BaseModel, Field, field_validator

# Configure Logging for safe API debugging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


# Target Gene Database with Market Research Metrics (Prevalence on a 1-10 scale)
MARKET_RESEARCH_GENES = {
    "PCSK9":   {"disease": "Hypercholesterolemia", "prevalence": 9.5},
    "HTT":     {"disease": "Huntington's Disease", "prevalence": 4.0},
    "HBB":     {"disease": "Sickle Cell / Thalassemia", "prevalence": 8.0},
    "TTR":     {"disease": "ATTR Amyloidosis", "prevalence": 5.5},
    "ANGPTL3": {"disease": "Dyslipidemia", "prevalence": 8.5},
    "DNMT1":   {"disease": "Hereditary Sensory Neuropathy", "prevalence": 2.5},
    "CFTR":    {"disease": "Cystic Fibrosis", "prevalence": 7.0}
}


class TrialRecord(BaseModel):
    """Pydantic model for strict field extraction and validation."""
    nct_id: str = Field(..., description="NCTId")
    brief_title: str = Field(..., description="BriefTitle")
    official_title: Optional[str] = "N/A"
    overall_status: str = Field(..., description="OverallStatus")
    phase: str = Field(..., description="Phase")
    lead_sponsor_name: str = Field(..., description="LeadSponsorName")
    conditions: List[str] = []
    interventions: List[str] = []

    @field_validator("nct_id", mode="before")
    def validate_nct(cls, v):
        return v if v and str(v).strip() != "" else "UNKNOWN_NCT"

    @field_validator("brief_title", "overall_status", "phase", "lead_sponsor_name", mode="before")
    def sanitize_strings(cls, v):
        if not v or str(v).strip() == "":
            return "Not Specified"
        return str(v).strip()


def extract_study_fields(study_data: dict) -> Optional[TrialRecord]:
    """Step 2: Safely extracts required NCTId, BriefTitle, OverallStatus, Phase, and LeadSponsorName."""
    try:
        protocol = study_data.get("protocolSection", {})
        
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
        design_module = protocol.get("designModule", {})
        conditions_module = protocol.get("conditionsModule", {})
        interventions_module = protocol.get("armsInterventionsModule", {})
        
        # Format Phase array to clean string
        phases = design_module.get("phases", ["Not Specified"])
        phase_str = ", ".join(phases) if isinstance(phases, list) else str(phases)
        
        # Interventions array extraction
        interventions_raw = interventions_module.get("interventions", [])
        intervention_names = [i.get("name", "") for i in interventions_raw if i.get("name")]

        return TrialRecord(
            nct_id=id_module.get("nctId"),
            brief_title=id_module.get("briefTitle"),
            official_title=id_module.get("officialTitle", "N/A"),
            overall_status=status_module.get("overallStatus"),
            phase=phase_str,
            lead_sponsor_name=sponsor_module.get("leadSponsor", {}).get("name"),
            conditions=conditions_module.get("conditions", []),
            interventions=intervention_names
        )
    except Exception as e:
        logging.warning(f"Error parsing study record: {e}")
        return None


def fetch_clinical_trials(keywords: List[str], max_pages_per_keyword: int = 2) -> pd.DataFrame:
    """Step 1: Robust API connection with timeouts and retry safety."""
    logging.info("[Step 1] Connecting to ClinicalTrials.gov API v2...")
    endpoint = "https://clinicaltrials.gov/api/v2/studies"
    records = []
    
    for term in keywords:
        page_token = None
        for page in range(max_pages_per_keyword):
            params = {"query.term": term, "pageSize": 50, "format": "json"}
            if page_token:
                params["pageToken"] = page_token
                
            try:
                response = requests.get(endpoint, params=params, timeout=12)
                response.raise_for_status()
                payload = response.json()
            except Exception as err:
                logging.error(f"Failed to fetch term '{term}' on page {page + 1}: {err}")
                break

            studies = payload.get("studies", [])
            for study in studies:
                record = extract_study_fields(study)
                if record:
                    records.append(record.model_dump())

            page_token = payload.get("nextPageToken")
            if not page_token:
                break

    df = pd.DataFrame(records)
    if not df.empty:
        df = df.drop_duplicates(subset=["nct_id"]).reset_index(drop=True)
        logging.info(f"[Step 2] Successfully extracted {len(df)} unique trial records.")
    else:
        logging.warning("No records fetched. Creating empty fallback DataFrame.")
        df = pd.DataFrame(columns=["nct_id", "brief_title", "official_title", "overall_status", "phase", "lead_sponsor_name", "conditions", "interventions"])
    
    return df



def fetch_patent_data(gene_targets: List[str]) -> Dict[str, int]:
    """Step 3: Queries Europe PMC REST API for patent literature saturation per gene."""
    logging.info("[Step 3] Gathering market research & patent saturation data...")
    patent_counts = {}
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    
    for gene in gene_targets:
        query_str = f'SRC:PAT AND ("CRISPR" OR "Gene Editing") AND "{gene}"'
        try:
            res = requests.get(url, params={"query": query_str, "format": "json", "pageSize": 1}, timeout=10)
            if res.status_code == 200:
                patent_counts[gene] = res.json().get("hitCount", 0)
            else:
                patent_counts[gene] = 0
        except Exception as e:
            logging.warning(f"Could not reach Europe PMC API for {gene}: {e}")
            patent_counts[gene] = 0

    return patent_counts


def analyze_market_opportunity(trials_df: pd.DataFrame, patent_counts: Dict[str, int]) -> pd.DataFrame:
    """Maps gene targets in trials and computes Opportunity Index."""
    def detect_genes(row):
        text = f"{row['brief_title']} {row['official_title']} {' '.join(row['conditions'])} {' '.join(row['interventions'])}"
        found = [g for g in MARKET_RESEARCH_GENES.keys() if re.search(r'\b' + re.escape(g) + r'\b', text, re.IGNORECASE)]
        return found if found else ["Uncategorized"]

    if not trials_df.empty:
        trials_df["detected_genes"] = trials_df.apply(detect_genes, axis=1)
        exploded = trials_df.explode("detected_genes")
        gene_trial_counts = exploded["detected_genes"].value_counts().to_dict()
    else:
        trials_df["detected_genes"] = [[]]
        gene_trial_counts = {}

    rows = []
    for gene, meta in MARKET_RESEARCH_GENES.items():
        active_trials = gene_trial_counts.get(gene, 0)
        patents = patent_counts.get(gene, 0)
        prev = meta["prevalence"]
        
        # Formula: Opportunity Index = Prevalence / ((Active Trials + 1) * log10(Patents + 10))
        opp_index = round(prev / ((active_trials + 1) * math.log10(patents + 10)), 2)
        
        rows.append({
            "Gene Target": gene,
            "Target Indication": meta["disease"],
            "Active Trials": active_trials,
            "Patent Count": patents,
            "Prevalence Index": prev,
            "Opportunity Score": opp_index
        })

    return pd.DataFrame(rows)



# (FOREST / BOTANICAL THEME)


def build_forest_dashboard(trials_df: pd.DataFrame, market_df: pd.DataFrame, output_filepath: str = "crispr_dashboard.html"):
    """Step 4: Renders a website formatted with custom HTML & CSS inspired by a forest motif."""
    logging.info("[Step 4] Building forest-themed HTML website...")

    forest_colors = ["#2A4735", "#52B788", "#74C69D", "#D4AF37", "#1B3B2B", "#85A389"]


    phase_counts = trials_df["phase"].value_counts().reset_index() if not trials_df.empty else pd.DataFrame({"phase": ["N/A"], "count": [0]})
    fig_phase = px.bar(
        phase_counts,
        x="phase",
        y="count",
        labels={"phase": "Trial Phase", "count": "Number of Studies"},
        title="Clinical Trial Pipeline Stage Distribution",
        color_discrete_sequence=["#52B788"]
    )
    fig_phase.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#E6EFE9", family="Georgia, serif"),
        margin=dict(l=20, r=20, t=50, b=20)
    )
    fig_phase.update_xaxes(showgrid=False, zeroline=False, color="#A3C1AD")
    fig_phase.update_yaxes(showgrid=True, gridcolor="#1B3B2B", color="#A3C1AD")

    fig_scatter = px.scatter(
        market_df,
        x="Patent Count",
        y="Active Trials",
        size="Prevalence Index",
        color="Opportunity Score",
        hover_name="Gene Target",
        text="Gene Target",
        title="IP Saturation vs. Clinical Activity",
        color_continuous_scale=["#1B3B2B", "#52B788", "#D4AF37"]
    )
    fig_scatter.update_traces(textposition="top center", marker=dict(line=dict(width=1, color="#D4AF37")))
    fig_scatter.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#E6EFE9", family="Georgia, serif"),
        margin=dict(l=20, r=20, t=50, b=20)
    )
    fig_scatter.update_xaxes(showgrid=True, gridcolor="#1B3B2B", color="#A3C1AD")
    fig_scatter.update_yaxes(showgrid=True, gridcolor="#1B3B2B", color="#A3C1AD")

    phase_chart_html = fig_phase.to_html(full_html=False, include_plotlyjs="cdn")
    scatter_chart_html = fig_scatter.to_html(full_html=False, include_plotlyjs="cdn")


    table_market = market_df.to_html(index=False, classes="forest-table")
    
    display_trials = trials_df[["nct_id", "brief_title", "overall_status", "phase", "lead_sponsor_name"]].head(12)
    display_trials.columns = ["NCTId", "BriefTitle", "OverallStatus", "Phase", "LeadSponsorName"]
    table_trials = display_trials.to_html(index=False, classes="forest-table")

    top_opportunity = market_df.loc[market_df["Opportunity Score"].idxmax()]["Gene Target"] if not market_df.empty else "N/A"

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Editing Market & Clinical Intelligence</title>
    <style>
        :root {{
            --bg-deep-forest: #0A1C14;
            --bg-card: #13271D;
            --border-moss: #224030;
            --text-main: #E8EFE9;
            --text-muted: #A3C1AD;
            --accent-gold: #D4AF37;
            --accent-emerald: #52B788;
        }}

        body {{
            background-color: var(--bg-deep-forest);
            color: var(--text-main);
            font-family: 'Segoe UI', Georgia, serif;
            margin: 0;
            padding: 30px;
            background-image: radial-gradient(circle at 50% 0%, #153324 0%, #0A1C14 75%);
        }}

        .container {{
            max-width: 1280px;
            margin: 0 auto;
        }}

        header {{
            border-bottom: 2px solid var(--border-moss);
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}

        h1 {{
            font-size: 2.2rem;
            color: var(--accent-gold);
            margin: 0 0 8px 0;
            font-weight: 400;
            letter-spacing: 1px;
        }}

        p.subtitle {{
            color: var(--text-muted);
            margin: 0;
            font-size: 1rem;
        }}

        .kpi-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}

        .kpi-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-moss);
            border-radius: 8px;
            padding: 20px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
        }}

        .kpi-value {{
            font-size: 2.2rem;
            font-weight: bold;
            color: var(--accent-emerald);
            margin-bottom: 4px;
        }}

        .kpi-label {{
            font-size: 0.85rem;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 1px;
        }}

        .content-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-moss);
            border-radius: 8px;
            padding: 24px;
            margin-bottom: 30px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
        }}

        h3 {{
            color: var(--accent-gold);
            margin-top: 0;
            font-weight: 400;
            border-bottom: 1px solid var(--border-moss);
            padding-bottom: 10px;
        }}

        .grid-2col {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }}

        @media (max-width: 900px) {{
            .grid-2col {{ grid-template-columns: 1fr; }}
        }}

        .forest-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            font-size: 0.95rem;
        }}

        .forest-table th {{
            background-color: #1A3628;
            color: var(--accent-gold);
            text-align: left;
            padding: 12px;
            border-bottom: 2px solid var(--border-moss);
        }}

        .forest-table td {{
            padding: 10px 12px;
            border-bottom: 1px solid var(--border-moss);
            color: var(--text-main);
        }}

        .forest-table tr:hover {{
            background-color: #193325;
        }}

        footer {{
            text-align: center;
            color: var(--text-muted);
            font-size: 0.85rem;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid var(--border-moss);
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Gene Editing Intelligence Dashboard</h1>
            <p class="subtitle">Clinical Trial Extraction (NCTId, Title, Status, Phase, Sponsor) & Market Opportunity Analysis</p>
        </header>

        <!-- KPI Summary Cards -->
        <div class="kpi-grid">
            <div class="kpi-card">
                <div class="kpi-value">{len(trials_df)}</div>
                <div class="kpi-label">Extracted Clinical Studies</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{len(MARKET_RESEARCH_GENES)}</div>
                <div class="kpi-label">Target Genes Tracked</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{market_df['Patent Count'].sum():,}</div>
                <div class="kpi-label">Indexed IP Patents</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{top_opportunity}</div>
                <div class="kpi-label">Highest Opportunity Target</div>
            </div>
        </div>

        <!-- Step 3: Market Research Table -->
        <div class="content-card">
            <h3>Target Gene Market Opportunity Index</h3>
            <p style="color: var(--text-muted); font-size: 0.9rem;">Calculated using disease prevalence relative to clinical saturation and patent density.</p>
            <div style="overflow-x: auto;">
                {table_market}
            </div>
        </div>

        <!-- Charts Grid -->
        <div class="grid-2col">
            <div class="content-card">
                {phase_chart_html}
            </div>
            <div class="content-card">
                {scatter_chart_html}
            </div>
        </div>

        <!-- Step 2: Extracted Trial Data Table -->
        <div class="content-card">
            <h3>Extracted Clinical Trial Data (ClinicalTrials.gov)</h3>
            <p style="color: var(--text-muted); font-size: 0.9rem;">Displaying extracted fields: NCTId, BriefTitle, OverallStatus, Phase, and LeadSponsorName.</p>
            <div style="overflow-x: auto;">
                {table_trials}
            </div>
        </div>

        <footer>
            Automated Gene Editing Pipeline &bull; ClinicalTrials.gov API & Europe PMC Patent API Integration
        </footer>
    </div>
</body>
</html>
"""

    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(html_template)

    logging.info(f"[Step 4 Complete] Website successfully generated: '{output_filepath}'")


# =====================================================================
# MAIN PIPELINE EXECUTION
# =====================================================================

if __name__ == "__main__":
    search_keywords = ["CRISPR", "Gene Editing", "Base Editing"]
    gene_targets = list(MARKET_RESEARCH_GENES.keys())

    # 1) API Connection & Search Working
    # 2) Data Extraction Working (NCTId, BriefTitle, OverallStatus, Phase, LeadSponsorName)
    trials_data = fetch_clinical_trials(search_keywords, max_pages_per_keyword=2)

    # 3) Market Research Gathered (Europe PMC Patents + Prevalence)
    patent_counts = fetch_patent_data(gene_targets)
    market_analytics = analyze_market_opportunity(trials_data, patent_counts)

    # Export structured raw datasets
    trials_data.to_csv("trials_dataset.csv", index=False)
    market_analytics.to_csv("opportunity_matrix.csv", index=False)

    # 4) Website Built Displaying the Data (Forest Aesthetic)
    build_forest_dashboard(trials_data, market_analytics, output_filepath="crispr_dashboard.html")
