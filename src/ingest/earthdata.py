import requests
from utils import json2xml, checkip
import json
import xmltodict
import logging

#

# global

class Earthdata(object):

    def __init__(self):
        self.token = None
        self.base_url = 'https://cmr.earthdata.nasa.gov'
        self.headers = {'Content-Type': 'application/xml', 'Client-Id': 'WRFXPY'}

    def login(self,username, password):
        # see https://wiki.earthdata.nasa.gov/display/CMR/CMR+Client+Partner+User+Guide#CMRClientPartnerUserGuide-CreatingaToken
        if self.token is not None:
            logging.warning('Already logged into %s' % self.base_url)
            return True
        url = self.base_url +'/legacy-services/rest/tokens'
        self.ip = checkip()
        data ={'token': {'username': username, 
            'password': password, 
            'user_ip_address': self.ip, 'client_id': 'WRFXPY'}}
        try:
            r = requests.post(url, data=json2xml(data), headers=headers)
        except Exception as e:
            logging.error(e)
            return False
        if r.status_code == 201:
            try:
                self.token = xmltodict.parse(r.text)['token']['id']
                logging.info('Successfully logged into %s, token=%s' % (self.base_url, self.token))
                return True
            except:
                logging.error('%s did not return a token' % self.base_url)
        else:
            logging.error('%s return code %s' % (self.base_url, r.status_code))
        logging.error('Cannot log into %s' % self.base_url) 
        return False
    
    def logout(self):
        if self.token is None:
            logging.warning('Already logged out of %s' % self.base_url)
        else:
            url = self.base_url +'/legacy-services/rest/tokens/' + self.token
            try:
                r = requests.delete(url, headers=headers)
                if r.status_code == 204:
                    self.token = None
                    logging.info('Successfully logged out of %s' % self.base_url)
                    return
            except Exception as e:
                logging.error(e)
            logging.warning('Failed to log out of %s' % self.base_url)

    # tentative stubs 

class MODIS(Earthdata):
    # special things for the MODIS instrument

class TERRA(MODIS):
    # anythnig special for the TERRA satellite

class AQUA(MODIS): 
    # anything special for the AQUA satellite

class VIIRS(Earthdata):
    # anything special for the VIIRS instrument

class SUOMI_NPP(VIIRS):
    # anything special for the SUOMI_NPP satellite

def search(): 
    # find all granules within a given box and time interval

def download_granules():
    # do thins

def extract_maps():
    # return list of rectangular arrays of pixel values with geolocation 
    # for each pixel: active fires detection status, confidence level, fire radiative power, longitude, latitude


