import requests, logging

class response_object(object):
    status_code = 0 
    def __init__(self,status_code):
        self.status_code = status_code


def readhead(url):
    try:
        ret=requests.head(url)
        ret.raise_for_status()
    except (requests.RequestException, requests.exceptions.Timeout, requests.exceptions.TooManyRedirects, requests.exceptions.ConnectionError, requests.exceptions.HTTPError, requests.exceptions.Timeout) as e:
        logging.warning(e)
        ret = response_object(-1)
    return ret 
