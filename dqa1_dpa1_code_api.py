import requests
import base64
import json

# Test script to send a request to UNOS Developer API to get DQA1 and DPA1 HLA field decoding
# HLA data is stored in numerical codes instead of IMGT/HLA string nomenclature

# VXM app
# https://developer.unos.org/my-apps/72b2149c-24cf-45e9-a4e0-e363d6941060


# credentials for app - save the key in hidden files later
tulane_app_id_file = open ("tulane_app_id.key")
tulane_app_id = ""
for line in tulane_app_id_file:
    tulane_app_id = line

beta_public_key_file = open ("unos_api_public.key")
beta_public_key = ""
for line in beta_public_key_file:
    beta_public_key = line

beta_secret_file = open ("unos_api_secret.key")
for line in beta_secret_file:
    beta_secret_key = line
client_credentials = beta_public_key + ":" + beta_secret_key

# base64 encode the client credentials and convert to raw - encode and decode ASCII
client_credentials_bytes = client_credentials.encode("ascii")
credentials_base64 = base64.b64encode(client_credentials_bytes).decode('ascii')

# create an authorization header as a dictionary
header_authorization = {'Authorization': credentials_base64}

# token request URL for DPA1 locus
# https://developer.unos.org/docs/histo-lab-management---test/1/routes/deceased-donor/v1/lookups/hla/dpa1-locus/get

token_request_URL = "https://api-beta.unos.org/oauth/accesstoken?grant_type=client_credentials"

# requests post to request access token
r = requests.post(token_request_URL, headers=header_authorization, data={})
print (r.content)

# extract access token from JSON response
token =  r.json()['access_token']
token_bearer = "Bearer " + token
print (token)

# use access token to request DQA1 and DPA1 URLs

DQA1_URL = "https://api-beta.unos.org/deceased-donor/v1/lookups/hla/dqa1-locus"
DPA1_URL = "https://api-beta.unos.org/deceased-donor/v1/lookups/hla/dpa1-locus"

header_authorization = {'Authorization': token_bearer }
header_multiple = {'Authorization': token_bearer,
                   'accept': 'application/json'}

response_bytes = requests.get(DPA1_URL, headers=header_multiple, data={})
print (response_bytes.content)

json_data = json.loads(response_bytes.content.decode('utf-8'))

formatted_json = json.dumps(json_data, indent=4)

print (formatted_json)

# b'{"fault":{"faultstring":"Invalid API call as no apiproduct match found","detail":{"errorcode":"keymanagement.service.InvalidAPICallAsNoApiProductMatchFound"}}}'