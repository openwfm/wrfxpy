import requests
from utils import json2xml, checkip
import json
import xmltodict
import logging

token = None

def login(username, password):
    # see https://wiki.earthdata.nasa.gov/display/CMR/CMR+Client+Partner+User+Guide#CMRClientPartnerUserGuide-CreatingaToken
    base_url = 'https://cmr.earthdata.nasa.gov'
    url = base_url +'/legacy-services/rest/tokens'
    data ={'token': {'username': username, 
        'password': password, 
        'user_ip_address': checkip(), 'client_id': 'WRFXPY'}}
    headers = {'Content-Type': 'application/xml'}
    try:
        r = requests.post(url, data=json2xml(data), headers=headers)
    except Exception as e:
        logging.error(e)
        return False
    if r.status_code == 201:
        try:
            token = xmltodict.parse(r.text)['token']['id']
            logging.info('Successfully logged into %s, token=%s' % (base_url, token))
            return True
        except:
            logging.error('%s did not return a token' % base_url)
    else:
        logging.error('%s return code %s' % (base_url, r.status_code))
    logging.error('Cannot log into %s' % base_url) 
    return False

