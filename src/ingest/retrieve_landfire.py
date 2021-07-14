from ingest.downloader import download_url
from geo.geo_utils import check_bbox
from utils import ensure_dir,move,remove,delete
import six.moves.urllib.request as request
import os.path as osp
import sys, time, zipfile, glob, random, logging

base_url = 'https://landfire.cr.usgs.gov/axis2/services/DownloadService'
fuel_opt = 'siz=2&key=FE0&ras=1&pfm=GeoTIFF&imsurl=-1&ms=-1&att=-1&lay=-1&fid=-1&dlpre=lf&lft={0}&rgt={1}&top={3}&bot={2}&wmd=1&mur=https://landfire.cr.usgs.gov/distmeta/servlet/gov.usgs.edc.MetaBuilder&mcd=FE0&mdf=HTML&arc=ZIP&sde=Landfire/US_200FBFM13&msd=LANDFIRE.US_140_SPAT_MASTER&zun=METERS&prj=102039&rsp=0&bnd=&bndnm=&csx=30.0&csy=30.0&ics=&ORIG=RVS'
elev_opt = 'siz=4&key=F0F&ras=1&pfm=GeoTIFF&imsurl=-1&ms=-1&att=-1&lay=-1&fid=-1&dlpre=lf&lft={0}&rgt={1}&top={3}&bot={2}&wmd=1&mur=https://landfire.cr.usgs.gov/distmeta/servlet/gov.usgs.edc.MetaBuilder&mcd=F0F&mdf=HTML&arc=ZIP&sde=Landfire/US_DEM2016&msd=LANDFIRE.US_140_SPAT_MASTER&zun=METERS&prj=102039&rsp=0&bnd=&bndnm=&csx=30.0&csy=30.0&ics=&ORIG=RVS'

class LandfireError(Exception):
    """
    Raised when a SatSource cannot retrieve satellite data.
    """
    pass

def server_service(service,args):
    url = osp.join(base_url,service+'?'+args)
    logging.info('server_service - service {0}, URL {1}'.format(service,url))
    sec = 5+random.random()*10
    time.sleep(sec)
    try:
        page = request.urlopen(url)
    except IOError as e:
        if hasattr(e, 'reason'):
            logging.error('failed to reach a server.')
            logging.error('reason: {}'.format(e.reason))
        elif hasattr(e, 'code'):
            logging.error('the server couldn\'t fulfill the request.')
            logging.error('code: ', e.code)
        raise LandfireError('Failed downloading the data...')
    else:
        result_page = str(page.read())
        if result_page.find('VALID>false') > -1:
            logging.error('problem with {0} request string {1}'.format(service,url))
            raise LandfireError('Failed downloading the data...')
        # remove carriage returns
        result_page = result_page.replace('&#xd;\\n',' ')
        # parse out return text
        startPos = result_page.find('<ns:return>') + 11
        endPos = result_page.find('</ns:return>')
        result = result_page[startPos:endPos]
    return result

def retrieve_landfire(opt_search,local_path,file_name):
    file_path = osp.join(local_path,file_name + '.tif')
    info_path = file_path + '.size'
    if osp.exists(file_path):
        if osp.getsize(file_path) == int(open(info_path).read()):
            logging.info('file {} already exists locally'.format(file_path))
            return file_path

    # submit a job to the Download Service
    requestID = server_service('initiateDownload',opt_search)
    logging.info('request ID: {}'.format(requestID))

    arg = 'downloadID={}'.format(requestID)
    stat = None
    retries = 50
    while stat != 400 and retries:
        # call Download service with request id to get status
        status = server_service('getDownloadStatus',arg)
        logging.info('status: {}'.format(status))
        stat = int(status.split(',')[0])
        retries -= 1

    if not retries and stat != 400:
        logging.error('maximum number of retries, the server is not responding...')
        raise LandfireError('Failed downloading the data...')

    # once a status of 400 has been received, retrieve from the URL
    url = server_service('getData',arg)
    logging.info('download URL: {}'.format(url))
    download_url(url,local_path+'.zip')

    # send complete message back to server so it can cleanup the job
    status = server_service('setDownloadComplete',arg)
    logging.info('status: {}'.format(status))

    # unzip and save
    local_zip = local_path+'.zip'
    local_tmp = osp.join(local_path,'tmp')
    with zipfile.ZipFile(local_zip,'r') as f:
        f.extractall(local_tmp)
    tif_files = glob.glob(osp.join(local_tmp,'*.tif'))
    if len(tif_files) == 1: 
        tif_file = tif_files[0]
        file_size = osp.getsize(tif_file)
        move(tif_file,file_path)
        with open(ensure_dir(info_path),'w') as f:
            f.write(str(file_size))
    else:
        logging.warning('not enought or too many TIF files, skipping...')
        return None
    remove(local_zip)
    remove(local_zip + '.size')
    delete(local_tmp)
    logging.info('data correctly retrieved as {}'.format(file_path))
    return file_path

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: ./retrieve_landfire.sh min_lon,max_lon,min_lat,max_lat')
        print('Example: ./retrieve_landfire.sh -112.8115,-112.1661,39.4820,39.9750')
        sys.exit(-1)
    else:
        coords = sys.argv[1].split(',')
        if len(coords) == 4:
            bbox = tuple([float(c) for c in coords])
        else:
            print('Invalid argument: {}'.format(sys.argv[1]))
            print('Usage: ./retrieve_landfire.sh min_lon,max_lon,min_lat,max_lat')
            print('Example: ./retrieve_landfire.sh -112.55,-112.4,39.65,39.8')
            sys.exit(-1)
        if not check_bbox(bbox):
            print('Invalid bounding box: {}'.format(sys.argv[1]))
            print('Usage: ./retrieve_landfire.sh min_lon,max_lon,min_lat,max_lat')
            print('Example: ./retrieve_landfire.sh -112.55,-112.4,39.65,39.8')
            sys.exit(-1)

    file_name = 'retrieve_landfire_{0:.4f}_{1:.4f}_{2:.4f}_{3:.4f}'.format(*bbox)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info('Retrieving LANDFIRE fuel data')
    retrieve_landfire(fuel_opt.format(*bbox),'ingest/landfire/fuel',file_name)
    logging.info('Retrieving LANDFIRE elevation data')
    retrieve_landfire(elev_opt.format(*bbox),'ingest/landfire/elevation',file_name)
