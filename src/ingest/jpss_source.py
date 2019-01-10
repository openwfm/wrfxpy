#
# Angel Farguell, CU Denver
#

import glob, re, datetime, logging, requests
import os.path as osp
import numpy as np
from urllib import urlopen
from cmr import GranuleQuery
from utils import Dict
from downloader import download_url, DownloadError

class JPSSError(Exception):
	"""
	Raised when a JPSSSource cannot retrieve JPSS data.
	"""
	pass

class JPSSSource(object):
	"""
	The parent class of all JPSS sources that implements common functionality, for example

	- local JPSS validation (file size check)
	- JPSS retrieval with retries (smart check whether server implements http-range)
	"""

	def __init__(self, js):
		"""
		Initialize jpss source with ingest directory (where JPSS files are stored).

		:param js: job structure with at least ingest_path root of JPSS storage and sys_install_path
		"""
		self.ingest_dir=osp.abspath(osp.join(js.get('ingest_path','ingest'),self.prefix))
		self.cache_dir=osp.abspath(js.get('cache_path','cache'))
		self.sys_dir=osp.abspath(js.get('sys_install_path',None))

	def available_locally_jpss(self, path):
		"""
		Check if a JPSS file is available locally and if it's file size checks out.

		:param path: the JPSS file path
		"""
		info_path = path + '.size'
		if osp.exists(path) and osp.exists(info_path):
			content_size = int(open(info_path).read())
			return osp.getsize(path) == content_size
		else:
			return False

	def search_api_jpss(self, prod, bbox, time, version=None):
		"""
		API search of the different satellite granules and return metadata dictionary

		:param prod: JPSS product
		:param bbox: polygon with the search bounding box
		:param time: time interval (init_time_iso,final_time_iso)

		:return metas: a dictionary with all the metadata for the API search
		"""
		maxg=100
		sname=self.prefix+prod
		api = GranuleQuery()
		if not version:
			search = api.parameters(
								short_name=sname,
								downloadable=True,
								polygon=bbox,
								temporal=time
								)
		else:
			search = api.parameters(
								short_name=sname,
								downloadable=True,
								polygon=bbox,
								temporal=time,
								version=version
								)
		sh=search.hits()
		logging.info("CMR API: %s gets %s hits in this range" % (sname,sh))
		if sh>maxg:
			logging.warning("The number of hits %s is larger than the limit %s." % (sh,maxg))
			logging.warning("Any JPSS data with prefix %s used." % self.prefix)
			logging.warning("Use a reduced bounding box or a reduced time interval.")
			metas = []
		else:
			metas = api.get(sh)
		return metas

	def search_archive_jpss(self, prod, time, geo_metas):
		"""
		Archive search of the different satellite granules

		:param prod: string of product with version, for instance: '5000/VNP09'
		:param time: time interval (init_time,final_time)
		:param geo_metas: granules of the geolocation metadata
		:return metas: dictionary with the metadata of the search
		"""
		ids=['.'.join(gm['producer_granule_id'].split('.')[1:3]) for gm in geo_metas] # satellite ids in bounding box
		metas=[]
		dts=[datetime.datetime.strptime(d,'%Y-%m-%dT%H:%M:%SZ') for d in time]
		delta=dts[1]-dts[0]
		nh=int(delta.total_seconds()/3600)
		dates=[dts[0]+datetime.timedelta(seconds=3600*k) for k in range(1,nh+1)]
		fold=np.unique(['%d/%03d' % (date.timetuple().tm_year,date.timetuple().tm_yday) for date in dates])
		urls=[self.url+'/'+self.prefix+prod+'/'+f for f in fold]
		for u in urls:
			js=requests.get(u+'.json').json()
			for j in js:
				arg=np.argwhere(np.array(ids)=='.'.join(j['name'].split('.')[1:3]))
				if arg.size:
					ar=arg[0][0]
					g=Dict(j)
					g.links=[{'href': u+'/'+g.name}]
					g.time_start=geo_metas[ar]['time_start']
					g.time_end=geo_metas[ar]['time_end']
					metas.append(g)
		logging.info('Archive: %s gets %s hits in this range' % (self.prefix+prod,len(metas)))
		return metas

	def get_metas_jpss(self, bbox, time):
		"""
		Get all the meta data for all the necessary products

		:param bbox: polygon with the search bounding box
		:param time: time interval (init_time_iso,final_time_iso)
		:param maxg: max number of granules to process
		:return granules: dictionary with the metadata of all the products
		"""
		metas=Dict([])
		metas.geo=self.search_api_jpss('03'+self.geo_extra_prefix,bbox,time,self.version)
		metas.fire=self.search_api_jpss('14',bbox,time)
		if not metas.fire:
			metas.fire=self.search_archive_jpss('14',time,metas.geo)
		metas.ref=self.search_api_jpss('09',bbox,time)
		if not metas.ref:
			metas.ref=self.search_archive_jpss('09',time,metas.geo)
		return metas

	def group_files_jpss(self):
		"""
		Agrupate the geolocation (03), fire (14) (, and reflectance (09)) files of a specific product in a path

		:return files: list of dictionaries with geolocation (03), fire (14) (, and reflectance (09)) file paths of the specific product
		"""
		path=self.ingest_dir
		reg=self.prefix

		files=[Dict({'geo':k}) for k in glob.glob(path+'/'+reg+'03*')]
		filesf=glob.glob(path+'/'+reg+'14*')
		filesr=glob.glob(path+'/'+reg+'09*')
		if len(filesf)>0:
			for f in filesf:
				mf=re.split("/",f)
				if mf is not None:
					m=mf[-1].split('.')
					if m is not None:
						for k,g in enumerate(files):
							mmf=re.split("/",g.geo)
							mm=mmf[-1].split('.')
							if mm[0][1]==m[0][1] and mm[1]+'.'+mm[2]==m[1]+'.'+m[2]:
								files[k].fire=f
		if len(filesr)>0:
			for f in filesr:
				mf=re.split("/",f)
				if mf is not None:
					m=mf[-1].split('.')
					if m is not None:
						for k,g in enumerate(files):
							mmf=re.split("/",g.geo)
							mm=mmf[-1].split('.')
							if mm[0][1]==m[0][1] and mm[1]+'.'+mm[2]==m[1]+'.'+m[2]:
								files[k].ref=f
		return files

	def download_jpss(self, url):
		"""
		Download a JPSS file from a JPSS service

		:param url: the URL of the file
		"""
		logging.info('downloading %s jpss data from %s' % (self.prefix, url))
		jpss_path = osp.join(self.ingest_dir,osp.basename(url))
		if self.available_locally_jpss(jpss_path):
			logging.info('%s is available locally' % jpss_path)
		else:
			try:
				download_url(url, jpss_path)
			except DownloadError as e:
				logging.error('%s cannot download jpss file %s' % (self.prefix, url))
				logging.warning('Pleae check %s for %s' % (self.info_url, self.info))
				raise JPSSError('JPSSSource: failed to download file %s' % url)

	def read_jpss(self, files, metas, bounds):
		"""
		Read all the JPSS files from a specific JPSS service

		:param file: dictionary with geolocation (03), fire (14) (, and reflectance (09)) file paths for JPSS service
		:param bounds: spatial bounds tuple (lonmin,lonmax,latmin,latmax)
		:return ret: dictionary with Latitude, Longitude and fire mask arrays read
		"""
		pass

	def read_files_jpss(self, file, bounds):
		"""
		Read a specific JPSS files from a specific JPSS service

		:param file: dictionary with geolocation (03), fire (14) (, and reflectance (09)) file paths for JPSS service
		:param bounds: spatial bounds tuple (lonmin,lonmax,latmin,latmax)
		:return ret: dictionary with Latitude, Longitude and fire mask arrays read
		"""
		pass

	def retrieve_data_jpss(self, bounds, time):
		"""
		Retrieve JPSS data in a bounding box coordinates and time interval

		:param bounds: polygon with the search bounding box
		:param time: time interval (init_time_iso,final_time_iso)
		:return data: dictonary with all the data
		"""
		lonmin,lonmax,latmin,latmax=bounds
		bbox=[(lonmin,latmax),(lonmin,latmin),(lonmax,latmin),(lonmax,latmax),(lonmin,latmax)]

		metas=self.get_metas_jpss(bbox,time)

		logging.info('retrieve_jpss: saved %s metas for %s JPSS service' % (len(metas),self.prefix))

		for k, mm in metas.items():
			logging.info('downloading %s products' % k)
			for m in mm:
				url=m['links'][0]['href']
				self.download_jpss(url)

		#group=self.group_files_jpss()
		#jpss_data=self.read_files_jpss()

		#return jpss_data

		return metas


	# instance variables
	info_url=None
	info=None
	url=None
	prefix=None
	geo_extra_prefix=''
	platform=None
	version=None

