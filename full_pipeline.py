import os
import re
import math
import logging
import webbrowser
import requests
import pandas as pd
import numpy as np
import plotly.express as px
from typing import List, Optional, Dict
from pydantic import BaseModel, Field, field_validator

# Configure Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# =====================================================================
# CONFIGURATION & MARKET DATA
# =====================================================================

MARKET_RESEARCH_GENES = {
    "PCSK9":   {"disease": "Hypercholesterolemia", "prevalence": 9.5},
    "HTT":     {"disease": "Huntington's Disease", "prevalence": 4.0},
    "HBB":     {"disease": "Sickle Cell / Thalassemia", "prevalence": 8.0},
    "TTR":     {"disease": "ATTR Amyloidosis", "prevalence": 5.5},
    "ANGPTL3": {"disease": "Dyslipidemia", "prevalence": 8.5},
    "DNMT1":   {"disease": "Hereditary Sensory Neuropathy", "prevalence": 2.5},
    "CFTR":    {"disease": "Cystic Fibrosis", "prevalence": 7.0}
}

# =====================================================================
# STEP 1 & 2: API CONNECTION & STRICT DATA EXTRACTION
# =====================================================================

class TrialRecord(BaseModel):
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
        return v if v and str(v).strip() != "" else "N/A"

    @field_validator("brief_title", "overall_status", "phase", "lead_sponsor_name", mode="before")
    def sanitize_strings(cls, v):
        if not v or str(v).strip() == "":
            return "Not Specified"
        return str(v).strip()


def extract_study_fields(study_data: dict) -> Optional[TrialRecord]:
    """Safely extracts NCTId, BriefTitle, OverallStatus, Phase, and LeadSponsorName."""
    try:
        protocol = study_data.get("protocolSection", {})
        
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        sponsor_module = protocol.get("sponsorCollaboratorsModule", {})
        design_module = protocol.get("designModule", {})
        conditions_module = protocol.get("conditionsModule", {})
        interventions_module = protocol.get("armsInterventionsModule", {})
        
        phases = design_module.get("phases", ["Not Specified"])
        phase_str = ", ".join(phases) if isinstance(phases, list) else str(phases)
        
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
        logging.warning(f"Parsing error: {e}")
        return None


def fetch_clinical_trials(keywords: List[str], max_pages_per_keyword: int = 2) -> pd.DataFrame:
    """Connects to ClinicalTrials.gov v2 API and extracts records."""
    logging.info("[Step 1/4] Querying ClinicalTrials.gov v2 API...")
    endpoint = "https://clinicaltrials.gov/api/v2/studies"
    records = []
    
    for term in keywords:
        page_token = None
        for page in range(max_pages_per_keyword):
            params = {
                "query.term": term, 
                "pageSize": 50, 
                "format": "json",
                "countTotal": "true"
            }
            if page_token:
                params["pageToken"] = page_token
                
            try:
                response = requests.get(endpoint, params=params, timeout=12)
                response.raise_for_status()
                payload = response.json()
            except Exception as err:
                logging.error(f"Failed to fetch term '{term}': {err}")
                break

            for study in payload.get("studies", []):
                record = extract_study_fields(study)
                if record:
                    records.append(record.model_dump())

            page_token = payload.get("nextPageToken")
            if not page_token:
                break

    df = pd.DataFrame(records)
    if not df.empty:
        df = df.drop_duplicates(subset=["nct_id"]).reset_index(drop=True)
        logging.info(f"[Step 2/4] Successfully extracted {len(df)} clinical study records.")
    else:
        df = pd.DataFrame(columns=["nct_id", "brief_title", "official_title", "overall_status", "phase", "lead_sponsor_name", "conditions", "interventions"])
    
    return df


# =====================================================================
# STEP 3: MARKET RESEARCH & PATENT DENSITY
# =====================================================================

def fetch_patent_data(gene_targets: List[str]) -> Dict[str, int]:
    """Queries Europe PMC REST API for patent literature saturation per gene."""
    logging.info("[Step 3/4] Fetching IP & patent counts via Europe PMC API...")
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
        except Exception:
            patent_counts[gene] = 0

    return patent_counts


def analyze_market_opportunity(trials_df: pd.DataFrame, patent_counts: Dict[str, int]) -> pd.DataFrame:
    """Calculates Opportunity Score per gene target."""
    def detect_genes(row):
        text = f"{row['brief_title']} {row['official_title']} {' '.join(row['conditions'])} {' '.join(row['interventions'])}"
        found = [g for g in MARKET_RESEARCH_GENES.keys() if re.search(r'\b' + re.escape(g) + r'\b', text, re.IGNORECASE)]
        return found if found else ["Uncategorized"]

    if not trials_df.empty:
        trials_df["detected_genes"] = trials_df.apply(detect_genes, axis=1)
        exploded = trials_df.explode("detected_genes")
        gene_trial_counts = exploded["detected_genes"].value_counts().to_dict()
    else:
        gene_trial_counts = {}

    rows = []
    for gene, meta in MARKET_RESEARCH_GENES.items():
        active_trials = gene_trial_counts.get(gene, 0)
        patents = patent_counts.get(gene, 0)
        prev = meta["prevalence"]
        
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


# =====================================================================
# STEP 4: HTML & JAVASCRIPT BOTANICAL DASHBOARD
# =====================================================================

def build_forest_dashboard(trials_df: pd.DataFrame, market_df: pd.DataFrame, output_filepath: str = "crispr_dashboard.html"):
    """Generates an interactive dashboard with botanical fonts and search filters."""
    logging.info("[Step 4/4] Generating HTML website dashboard...")

    # Plotly Chart 1: Pipeline Stage Distribution
    phase_counts = trials_df["phase"].value_counts().reset_index() if not trials_df.empty else pd.DataFrame({"phase": ["N/A"], "count": [0]})
    fig_phase = px.bar(
        phase_counts,
        x="phase",
        y="count",
        labels={"phase": "Trial Phase", "count": "Studies"},
        title="Clinical Trial Pipeline Stage Distribution",
        color_discrete_sequence=["#52B788"]
    )
    fig_phase.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#E8EFE9", family="'Plus Jakarta Sans', sans-serif"),
        title_font=dict(family="'Cinzel', serif", color="#D4AF37", size=18),
        margin=dict(l=20, r=20, t=50, b=20)
    )
    fig_phase.update_xaxes(showgrid=False, color="#A3C1AD")
    fig_phase.update_yaxes(showgrid=True, gridcolor="#1B3B2B", color="#A3C1AD")

    # Plotly Chart 2: Patent Density vs Clinical Activity
    fig_scatter = px.scatter(
        market_df,
        x="Patent Count",
        y="Active Trials",
        size="Prevalence Index",
        color="Opportunity Score",
        hover_name="Gene Target",
        text="Gene Target",
        title="IP Saturation vs. Clinical Trial Density",
        color_continuous_scale=["#1B3B2B", "#52B788", "#D4AF37"]
    )
    fig_scatter.update_traces(textposition="top center", marker=dict(line=dict(width=1, color="#D4AF37")))
    fig_scatter.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#E8EFE9", family="'Plus Jakarta Sans', sans-serif"),
        title_font=dict(family="'Cinzel', serif", color="#D4AF37", size=18),
        margin=dict(l=20, r=20, t=50, b=20)
    )
    fig_scatter.update_xaxes(showgrid=True, gridcolor="#1B3B2B", color="#A3C1AD")
    fig_scatter.update_yaxes(showgrid=True, gridcolor="#1B3B2B", color="#A3C1AD")

    phase_chart_html = fig_phase.to_html(full_html=False, include_plotlyjs="cdn")
    scatter_chart_html = fig_scatter.to_html(full_html=False, include_plotlyjs="cdn")

    # Prepare Market Table HTML
    table_market_rows = ""
    for _, row in market_df.iterrows():
        table_market_rows += f"""
        <tr>
            <td><strong>{row['Gene Target']}</strong></td>
            <td>{row['Target Indication']}</td>
            <td>{row['Active Trials']}</td>
            <td>{row['Patent Count']:,}</td>
            <td>{row['Prevalence Index']}</td>
            <td style="color: #D4AF37; font-weight: bold;">{row['Opportunity Score']}</td>
        </tr>
        """

    # Prepare Trial Rows with JavaScript Search Support
    display_trials = trials_df[["nct_id", "brief_title", "overall_status", "phase", "lead_sponsor_name"]]
    trial_rows_html = ""
    for _, row in display_trials.iterrows():
        status_color = "#52B788" if "RECRUITING" in str(row['overall_status']).upper() else "#A3C1AD"
        trial_rows_html += f"""
        <tr class="trial-row">
            <td><code style="color: #74C69D;">{row['nct_id']}</code></td>
            <td>{row['brief_title']}</td>
            <td><span class="badge" style="border: 1px solid {status_color}; color: {status_color};">{row['overall_status']}</span></td>
            <td>{row['phase']}</td>
            <td>{row['lead_sponsor_name']}</td>
        </tr>
        """

    top_opportunity = market_df.loc[market_df["Opportunity Score"].idxmax()]["Gene Target"] if not market_df.empty else "N/A"

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Editing Intelligence & Market Dashboard</title>
    <!-- Google Fonts Integration -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Cinzel:wght@500;700&family=Plus+Jakarta+Sans:wght@300;400;600;700&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-deep-forest: #0B1912;
            --bg-card: #13271D;
            --border-moss: #224030;
            --text-main: #E8EFE9;
            --text-muted: #A3C1AD;
            --accent-gold: #D4AF37;
            --accent-emerald: #52B788;
        }}

        * {{
            box-sizing: border-box;
        }}

        body {{
            background-color: var(--bg-deep-forest);
            color: var(--text-main);
            font-family: 'Plus Jakarta Sans', sans-serif;
            margin: 0;
            padding: 30px;
            background-image: radial-gradient(circle at 50% 0%, #173827 0%, #0B1912 80%);
        }}

        .container {{
            max-width: 1300px;
            margin: 0 auto;
        }}

        header {{
            border-bottom: 1px solid var(--border-moss);
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}

        h1 {{
            font-family: 'Cinzel', serif;
            font-size: 2.4rem;
            color: var(--accent-gold);
            margin: 0 0 6px 0;
            letter-spacing: 1px;
        }}

        .subtitle {{
            color: var(--text-muted);
            margin: 0;
            font-size: 0.95rem;
            font-weight: 300;
        }}

        .kpi-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}

        .kpi-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-moss);
            border-radius: 10px;
            padding: 20px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.25);
        }}

        .kpi-value {{
            font-size: 2.2rem;
            font-weight: 700;
            color: var(--accent-emerald);
            margin-bottom: 4px;
        }}

        .kpi-label {{
            font-size: 0.8rem;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 1px;
        }}

        .content-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-moss);
            border-radius: 10px;
            padding: 24px;
            margin-bottom: 30px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.25);
        }}

        h2, h3 {{
            font-family: 'Cinzel', serif;
            color: var(--accent-gold);
            margin-top: 0;
            font-weight: 600;
        }}

        .grid-2col {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }}

        @media (max-width: 900px) {{
            .grid-2col {{ grid-template-columns: 1fr; }}
        }}

        /* Table Styling */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            font-size: 0.9rem;
        }}

        th {{
            font-family: 'Cinzel', serif;
            background-color: #1A3628;
            color: var(--accent-gold);
            text-align: left;
            padding: 12px;
            border-bottom: 2px solid var(--border-moss);
        }}

        td {{
            padding: 12px;
            border-bottom: 1px solid var(--border-moss);
            color: var(--text-main);
        }}

        tr:hover {{
            background-color: #183325;
        }}

        .badge {{
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 0.75rem;
            font-weight: 600;
        }}

        /* Interactive Controls */
        .search-box {{
            width: 100%;
            padding: 12px 16px;
            background: #0B1912;
            border: 1px solid var(--border-moss);
            border-radius: 6px;
            color: var(--text-main);
            font-family: 'Plus Jakarta Sans', sans-serif;
            font-size: 0.95rem;
            margin-bottom: 15px;
            outline: none;
        }}

        .search-box:focus {{
            border-color: var(--accent-emerald);
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
            <h1>Gene Editing Market & Clinical Intelligence</h1>
            <p class="subtitle">Extracting NCTId, BriefTitle, OverallStatus, Phase, and LeadSponsorName via ClinicalTrials.gov API</p>
        </header>

        <!-- KPI Summary Bar -->
        <div class="kpi-grid">
            <div class="kpi-card">
                <div class="kpi-value">{len(trials_df)}</div>
                <div class="kpi-label">Studies Extracted</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{len(MARKET_RESEARCH_GENES)}</div>
                <div class="kpi-label">Gene Targets Tracked</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{market_df['Patent Count'].sum():,}</div>
                <div class="kpi-label">Total IP Patents</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-value">{top_opportunity}</div>
                <div class="kpi-label">Top Market Target</div>
            </div>
        </div>

        <!-- Step 3: Market Opportunity Matrix -->
        <div class="content-card">
            <h3>Gene Target Opportunity Index</h3>
            <p style="color: var(--text-muted); font-size: 0.85rem; margin-bottom: 15px;">
                Opportunity Index = Disease Prevalence Score / [(Active Trials + 1) × log10(Patents + 10)]
            </p>
            <div style="overflow-x: auto;">
                <table>
                    <thead>
                        <tr>
                            <th>Gene Target</th>
                            <th>Target Indication</th>
                            <th>Active Trials</th>
                            <th>Patent Count</th>
                            <th>Prevalence Index</th>
                            <th>Opportunity Score</th>
                        </tr>
                    </thead>
                    <tbody>
                        {table_market_rows}
                    </tbody>
                </table>
            </div>
        </div>

        <!-- Visual Analytics -->
        <div class="grid-2col">
            <div class="content-card">
                {phase_chart_html}
            </div>
            <div class="content-card">
                {scatter_chart_html}
            </div>
        </div>

        <!-- Step 2: Extracted Clinical Trial Data with Live Filter -->
        <div class="content-card">
            <h3>Extracted Clinical Trials</h3>
            <input type="text" id="trialSearch" class="search-box" placeholder="Type to filter trials by keyword, status, sponsor, or NCT ID..." onkeyup="filterTrials()">
            <div style="overflow-x: auto;">
                <table id="trialsTable">
                    <thead>
                        <tr>
                            <th>NCTId</th>
                            <th>BriefTitle</th>
                            <th>OverallStatus</th>
                            <th>Phase</th>
                            <th>LeadSponsorName</th>
                        </tr>
                    </thead>
                    <tbody>
                        {trial_rows_html}
                    </tbody>
                </table>
            </div>
        </div>

        <footer>
            Botanical Market Dashboard &bull; ClinicalTrials.gov API v2 & Europe PMC Integration
        </footer>
    </div>

    <!-- JavaScript Live Search -->
    <script>
        function filterTrials() {{
            const query = document.getElementById('trialSearch').value.toLowerCase();
            const rows = document.querySelectorAll('#trialsTable tbody tr');
            
            rows.forEach(row => {{
                const text = row.innerText.toLowerCase();
                row.style.display = text.includes(query) ? '' : 'none';
            }});
        }}
    </script>
</body>
</html>
"""

    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(html_template)

    logging.info(f"[Step 4 Complete] Dashboard saved to '{output_filepath}'.")
    
    # Automatically open in browser
    try:
        webbrowser.open("file://" + os.path.realpath(output_filepath))
    except Exception:
        pass


# =====================================================================
# MAIN PIPELINE EXECUTION
# =====================================================================

if __name__ == "__main__":
    search_keywords = ["CRISPR", "Gene Editing", "Base Editing"]
    gene_targets = list(MARKET_RESEARCH_GENES.keys())

    # Step 1 & 2: API Search & Field Extraction
    trials_df = fetch_clinical_trials(search_keywords, max_pages_per_keyword=2)

    # Step 3: Market Research
    patent_counts = fetch_patent_data(gene_targets)
    market_df = analyze_market_opportunity(trials_df, patent_counts)

    # Export Datasets
    trials_df.to_csv("trials_dataset.csv", index=False)
    market_df.to_csv("opportunity_matrix.csv", index=False)

    # Step 4: Build & Launch Forest/Botanical Dashboard
    build_forest_dashboard(trials_df, market_df, output_filepath="crispr_dashboard.html")