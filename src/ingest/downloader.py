# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# requests library doesn't support ftp downloads
# they do have requests-ftp, but they do not have confidence in it
# "This library was cowboyed together in about 4 hours of total work,
#  has no tests, and relies on a few ugly hacks."
import requests
import urllib2

import os.path as osp
import logging

from utils import ensure_dir

class DownloadError(Exception):
    """
    Raised when the downloader is unable to retrieve a URL.
    """
    pass

def get_dList(url):
    dList = urllib2.urlopen(url).read().splitlines()
    listing = []
    for l in dList:
        listing.append(l.split()[-1])
    return listing

def download_url(url, local_path, max_retries=3):
    """
    Download a remote URL to the location local_path with retries.

    On download, the file size is first obtained and stored.  When the download completes,
    the file size is compared to the stored file.  This prevents broken downloads from
    contaminating the processing chain.

    :param url: the remote URL
    :param local_path: the path to the local file
    :param max_retries: how many times we may retry to download the file
    """

    logging.info('download_url %s as %s' % (url,local_path))

    if url[:6] == 'ftp://':
        # do ftp download to the local path
        r = urllib2.urlopen(url)

        content_size = int(r.headers['Content-Length'])
        info_path = local_path + '.size'
        open(ensure_dir(info_path), 'w').write(str(content_size))

        current_retries = max_retries

        while current_retries > 0:
            with open(ensure_dir(local_path),'wb') as f:
                while True:
                    chunk = r.read(1024*1024)
                    if not chunk:
                        break
                    f.write(chunk)

            if content_size != osp.getsize(local_path):
                os.remove(local_path)
                logging.info('download_url trying again, retries available %d' % retries_available)
                if current_retries == 0:
                    os.remove(local_path)
                    os.remove(info_path)
                    raise DownloadError('failed to download file %s' % url)
                current_retries -= 1

            else:
                break


    else:
        r = requests.get(url, stream=True)
        content_size = int(r.headers['Content-Length'])

        # dump the correct file size to an info file next to the grib file
        # when re-using the GRIB2 file, we check its file size against this record
        # to avoid using partial files
        info_path = local_path + '.size'
        open(ensure_dir(info_path), 'w').write(str(content_size))

        # stream the download to file
        with open(ensure_dir(local_path), 'wb') as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)

                # does the server accept byte range queries? e.g. the NOMADs server does
                accepts_ranges = 'bytes' in r.headers.get('Accept-Ranges', '')

                retries_available = max_retries
                file_size = osp.getsize(local_path)
                while file_size < content_size:
                    if retries_available > 0:
                        logging.info('download_url trying again, retries available %d' % retries_available)
                        if accepts_ranges:
                            # if range queries are supported, try to download only the missing portion of the file
                            headers = { 'Range' : 'bytes=%d-%d' % (file_size, content_size) }
                            r = requests.get(url, headers=headers, stream=True)
                            with open(local_path, 'ab') as f:
                                for chunk in r.iter_content(1024 * 1024):
                                    f.write(chunk)
                                    retries_available -= 1
                                    file_size = osp.getsize(local_path)
                                else:
                                    # call the entire function recursively, this will attempt to redownload the entire file
                                    # and overwrite previously downloaded data
                                    self.download_grib(url, local_path, max_retries-1)
                        else:
                            os.remove(local_path)
                            os.remove(info_path)
                            raise DownloadError('failed to download file %s' % url)


