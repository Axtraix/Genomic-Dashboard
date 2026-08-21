# Import necessary modules: Flask for web serving, requests for fetching API data
from flask import Flask, render_template
import requests

# Initialize the Flask web application
app = Flask(__name__)

# Define the root route ("/") that triggers when visiting the website
@app.route("/")
def dashboard():
   # Base API endpoint for fetching clinical study data
    url = "https://clinicaltrials.gov/api/v2/studies"
    # Query parameters to filter results (recruiting diabetes studies, limit to 5)
    params = {
        "query.cond": "Diabetes",
        "filter.overallStatus": "RECRUITING",
        "pageSize": 5
    }

    # Send a GET request to the ClinicalTrials.gov API
    response = requests.get(url, params=params)
    # List to hold processed study details for the frontend
    studies = []

    # Check if the API request was successful (HTTP status 200)
    if response.status_code == 200:
        data = response.json()
        raw_studies = data.get('studies', [])
        # Parse the JSON response body
        for s in raw_studies:
            protocol = s.get('protocolSection', {})
            id_mod = protocol.get('identificationModule', {})
            status_mod = protocol.get('statusModule', {})

            # Extract relevant fields and append them as a dictionary to our list
            studies.append({
                "nct_id": id_mod.get('nctId', 'N/A'),
                "title": id_mod.get('briefTitle', 'N/A'),
                "status": status_mod.get('overallStatus', 'N/A')
            })

# Pass the processed studies list to the HTML template for rendering
    return render_template("dashboard.html", studies=studies)

if __name__ == "__main__":
    app.run(debug=True)