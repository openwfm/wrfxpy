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

import requests
import os.path as osp

from utils import ensure_dir


class DownloadError(Exception):
    """
    Raised when the downloader is unable to retrieve a URL.
    """
    pass


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

