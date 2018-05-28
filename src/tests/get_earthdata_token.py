import requests
from utils import json2xml
import json
import xmltodict


def get_earthdata_token(username, password):
    url = 'https://cmr.earthdata.nasa.gov/legacy-services/rest/tokens'
    data ={'token': {'username': username, 
        'password': password, 
        'user_ip_address': '127.0.0.1', 'client_id': '12345'}}
    headers = {'Content-Type': 'application/xml'}
    r = requests.post(url, data=json2xml(data), headers=headers)
    return  xmltodict.parse(r.text)['token']['id']

