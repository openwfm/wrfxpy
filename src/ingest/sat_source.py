#
# Angel Farguell, CU Denver
#

import glob, re, datetime, logging, requests
import os.path as osp
import numpy as np
from urllib import urlopen
from cmr import GranuleQuery
from utils import Dict, utc_to_utcf
from downloader import download_url, DownloadError

class SatError(Exception):
	"""
	Raised when a SatSource cannot retrieve satellite data.
	"""
	pass

class SatSource(object):
	"""
	The parent class of all satellite sources that implements common functionality, for example

	- local satellite validation (file size check)
	- Satellite retrieval with retries (smart check whether server implements http-range)
	"""

	def __init__(self, js):
		"""
		Initialize satellite source with ingest directory (where satellite files are stored).

		:param js: job structure with at least ingest_path root of satellite storage and sys_install_path
		"""
		self.ingest_dir=osp.abspath(osp.join(js.get('ingest_path','ingest'),self.prefix))
		self.cache_dir=osp.abspath(js.get('cache_path','cache'))
		self.sys_dir=osp.abspath(js.get('sys_install_path',None))

	def available_locally_sat(self, path):
		"""
		Check if a satellite file is available locally and if it's file size checks out.

		:param path: the satellite file path
		"""
		info_path = path + '.size'
		if osp.exists(path) and osp.exists(info_path):
			content_size = int(open(info_path).read())
			return osp.getsize(path) == content_size
		else:
			return False

	def search_api_sat(self, sname, bbox, time, version=None):
		"""
		API search of the different satellite granules and return metadata dictionary

		:param sname: short name satellite product, ex: 'MOD03'
		:param bbox: polygon with the search bounding box
		:param time: time interval (init_time_iso,final_time_iso)

		:return metas: a dictionary with all the metadata for the API search
		"""
		maxg=100
		time_utcf=(utc_to_utcf(time[0]),utc_to_utcf(time[1]))
		api = GranuleQuery()
		if not version:
			search = api.parameters(
								short_name=sname,
								downloadable=True,
								polygon=bbox,
								temporal=time_utcf
								)
		else:
			search = api.parameters(
								short_name=sname,
								downloadable=True,
								polygon=bbox,
								temporal=time_utcf,
								version=version
								)
		sh=search.hits()
		logging.info("CMR API: %s gets %s hits in this range" % (sname,sh))
		if sh>maxg:
			logging.warning("The number of hits %s is larger than the limit %s." % (sh,maxg))
			logging.warning("Any satellite data with prefix %s used." % self.prefix)
			logging.warning("Use a reduced bounding box or a reduced time interval.")
			metas = []
		else:
			metas = api.get(sh)
		return metas

	def search_archive_sat(self, prod, time, geo_metas):
		"""
		Archive search of the different satellite granules

		:param prod: string of short name product, ex: 'VNP09'
		:param time: time interval (init_time_datetime,final_time_datetime)
		:param geo_metas: granules of the geolocation metadata
		:return metas: dictionary with the metadata of the search
		"""
		ids=['.'.join(gm['producer_granule_id'].split('.')[1:3]) for gm in geo_metas] # satellite ids in bounding box
		metas=[]
		delta=time[1]-time[0]
		nh=int(delta.total_seconds()/3600)
		dates=[time[0]+datetime.timedelta(seconds=3600*k) for k in range(1,nh+1)]
		fold=np.unique(['%d/%03d' % (date.timetuple().tm_year,date.timetuple().tm_yday) for date in dates])
		urls=[self.url+'/'+prod+'/'+f for f in fold]
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
					g.producer_granule_id=j['name']
					metas.append(g)
		logging.info('Archive: %s gets %s hits in this range' % (self.prefix+prod,len(metas)))
		return metas

	def get_metas_sat(self, bbox, time):
		"""
		Get all the meta data for all the necessary products

		:param bbox: polygon with the search bounding box
		:param time: time interval (init_time_datetime,final_time_datetime)
		:param maxg: max number of granules to process
		:return granules: dictionary with the metadata of all the products
		"""
		metas=Dict([])
		metas.geo=self.search_api_sat(self.geo_prefix,bbox,time,self.version)
		if not metas.geo:
			if not self.pre_geo_prefix:
				logging.warning('any geolocation data matches the search.')
				return metas
			else:
				pre_geo=self.search_api_sat(self.pre_geo_prefix,bbox,time,self.version)
				metas.geo=self.search_archive_sat(self.geo_prefix,time,pre_geo)
		metas.fire=self.search_api_sat(self.fire_prefix,bbox,time)
		if not metas.fire:
			metas.fire=self.search_archive_sat(self.fire_prefix,time,metas.geo)
		metas.ref=self.search_api_sat(self.ref_prefix,bbox,time)
		if not metas.ref:
			metas.ref=self.search_archive_sat(self.ref_prefix,time,metas.geo)
		return metas

	def download_sat(self, url):
		"""
		Download a satellite file from a satellite service

		:param url: the URL of the file
		"""
		logging.info('downloading %s satellite data from %s' % (self.prefix, url))
		sat_name = osp.basename(url)
		sat_path = osp.join(self.ingest_dir,sat_name)
		if self.available_locally_sat(sat_path):
			logging.info('%s is available locally' % sat_path)
			return {'url': url,'local_path': sat_path}
		else:
			try:
				download_url(url, sat_path)
				return {'url': url,'local_path': sat_path,'downloaded': datetime.datetime.now}
			except DownloadError as e:
				logging.error('%s cannot download satellite file %s' % (self.prefix, url))
				logging.warning('Pleae check %s for %s' % (self.info_url, self.info))
				return {}
				raise SatError('SatSource: failed to download file %s' % url)

	def retrieve_data_sat(self, bounds, from_utc, to_utc):
		"""
		Retrieve satellite data in a bounding box coordinates and time interval

		:param bounds: polygon with the search bounding box
		:param time: time interval (init_time_iso,final_time_iso)
		:return data: dictonary with all the data
		"""
		lonmin,lonmax,latmin,latmax = bounds
		bbox = [(lonmin,latmax),(lonmin,latmin),(lonmax,latmin),(lonmax,latmax),(lonmin,latmax)]
		time = (from_utc, to_utc)

		metas = self.get_metas_sat(bbox,time)

		logging.info('retrieve_data_sat: saved %s metas for %s satellite service' % (sum([len(m) for m in metas]),self.prefix))

		manifest = Dict({})

		for k, mm in metas.items():
			logging.info('downloading %s products' % k)
			for m in mm:
				id = m['producer_granule_id']
				url = m['links'][0]['href']
				m.update(self.download_sat(url))
				manifest.update({id: m})

		return manifest

	def group_files_sat(self):
		"""
		Agrupate the geolocation (03), fire (14) (, and reflectance (09)) files of a specific product in a path

		:return files: list of dictionaries with geolocation (03), fire (14) (, and reflectance (09)) file paths of the specific product
		"""
		path=self.ingest_dir

		files=[Dict({'geo':k}) for k in glob.glob(path+'/'+self.geo_prefix)]
		filesf=glob.glob(path+'/'+self.fire_prefix)
		filesr=glob.glob(path+'/'+self.ref_prefix)
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

	def read_sat(self, files, metas, bounds):
		"""
		Read all the satellite files from a specific satellite service

		:param file: dictionary with geolocation (03), fire (14) (, and reflectance (09)) file paths for satellite service
		:param bounds: spatial bounds tuple (lonmin,lonmax,latmin,latmax)
		:return ret: dictionary with Latitude, Longitude and fire mask arrays read
		"""
		pass

	def read_files_sat(self, file, bounds):
		"""
		Read a specific satellite files from a specific satellite service

		:param file: dictionary with geolocation (03), fire (14) (, and reflectance (09)) file paths for satellite service
		:param bounds: spatial bounds tuple (lonmin,lonmax,latmin,latmax)
		:return ret: dictionary with Latitude, Longitude and fire mask arrays read
		"""
		pass


	# instance variables
	info_url=None
	info=None
	url=None
	prefix=None
	pre_geo_prefix=None
	geo_prefix=None
	fire_prefix=None
	ref_prefix=None
	platform=None
	version=None

