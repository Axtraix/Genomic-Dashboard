import requests


# Base endpoint/ knowing the api
url= "https://clinicaltrails.gov/api/v2/studies"

#Define query params
params= {
    "query.cond": "Diabetes",
    "filter.overallStatus": "RECRUITING",
    "pageSize": 5
}

# Send GET request
response = requests.get(url,params=params)

#Check for success 
if response.status_code == 200:
    data = response.json()
    print(f"Total Studies Returned: {len(data.get('studies',[]))}")
else:
    print(f"Failed with status code: {response.status_code}")