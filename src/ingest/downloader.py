# Copyright (C) 2013-2016 Martin Vejmelka, CU Denver
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

import six.moves.urllib.request, six.moves.urllib.error, six.moves.urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

import os.path as osp
import os
import logging
import time
import subprocess
import random

from utils import ensure_dir, load_sys_cfg, remove

cfg = load_sys_cfg()
sleep_seconds_def = cfg.get('sleep_seconds', 20)
max_retries_def = cfg.get('max_retries', 3)
wget = cfg.get('wget','/usr/bin/wget')
wget_options=cfg.get('wget_options',["--read-timeout=1"])
download_sleep_seconds=cfg.get('download_sleep_seconds', 5)


class DownloadError(Exception):
    """
    Raised when the downloader is unable to retrieve a URL.
    """
    pass


def get_dList(url):
    dList = six.moves.urllib.request.urlopen(url).read().splitlines()
    listing = []
    for l in dList:
        listing.append(l.split()[-1])
    return listing


def request_url(url, token):
    use_urllib2 = url[:6] == 'ftp://'
    use_aws = url[:5] == 's3://'
    if token:
        if use_urllib2:
            r = six.moves.urllib.request.urlopen(six.moves.urllib.request.Request(url,headers={"Authorization": "Bearer %s" % token}))
        elif use_aws:
            r = subprocess.check_output('aws s3 ls {}'.format(url), shell=True)
        else:
            r = requests.get(url, stream=True, headers={"Authorization": "Bearer %s" % token})
    else:
        if use_urllib2:
            r = six.moves.urllib.request.urlopen(url)
        elif use_aws:
            r = subprocess.check_output('aws s3 ls {} --no-sign-request'.format(url), shell=True)
        else:
            r = requests.get(url, stream=True)
    if use_aws:
        lines = r.decode().split('\n')
        for line in lines:
            if line.endswith(osp.basename(url)):
                content_size = int(line.split()[2])
                if content_size == 0:
                    raise DownloadError('download_url content size is equal to 0')
    else:
        if 'content-length' in r.headers.keys():
            content_size = int(r.headers.get('content-length',0))
            if content_size == 0:
                raise DownloadError('download_url content size is equal to 0')
        else:
            logging.warning('download_url not header contet-length information')
    return content_size


def download_url(url, local_path, max_retries=max_retries_def, sleep_seconds=sleep_seconds_def, token=None, min_size=1):
    """
    Download a remote URL to the location local_path with retries.

    On download, the file size is first obtained and stored.  When the download completes,
    the file size is compared to the stored file.  This prevents broken downloads from
    contaminating the processing chain.

    :param url: the remote URL
    :param local_path: the path to the local file
    :param max_retries: how many times we may retry to download the file
    :param token: key to use for the download or None if not
    """
    logging.info('download_url - %s as %s' % (url, local_path))
    logging.debug('download_url - if download fails, will try %d times and wait %d seconds each time' % (max_retries, sleep_seconds))
    sec = random.random() * download_sleep_seconds
    logging.info('download_url -  sleeping %s seconds' % sec)
    time.sleep(sec)

    use_aws = url[:5] == 's3://'

    try:
        content_size = request_url(url, token)
    except Exception as e:
        if max_retries > 0:
            logging.info('donwload_url - not found, download_url trying again, retries available %d' % max_retries)
            logging.warning(e)
            logging.info('download_url - sleeping %s seconds' % sleep_seconds)
            time.sleep(sleep_seconds)
            download_url(url, local_path, max_retries = max_retries-1, token = token, min_size = min_size)
        return

    logging.info('download_url %s as %s' % (url,local_path))
    remove(local_path)
    if use_aws:
        command=['aws','s3','cp',url,ensure_dir(local_path)]
        if token is None:
            command.append('--no-sign-request')
    else:
        command=[wget,'-O',ensure_dir(local_path),url]
        for opt in wget_options:
            command.insert(1,opt)
        command.insert(1,'--tries={}'.format(max_retries_def))
        if token:
            command.insert(1,'--header=\'Authorization: Bearer %s\'' % token)
    logging.info(' '.join(command))
    subprocess.call(' '.join(command),shell=True)

    file_size = osp.getsize(local_path)

    # content size may have changed during download
    content_size = request_url(url, token)
    logging.info('download_url - local file size {} remote content size {} minimum size {}'.format(file_size, content_size, min_size))

    # it should be != but for some reason content_size is wrong sometimes
    if content_size is not None and int(file_size) < int(content_size) or int(file_size) < min_size:
        logging.warning('download_url - wrong file size, download_url trying again, retries available %d' % max_retries)
        if max_retries > 0:
            # call the entire function recursively, this will attempt to redownload the entire file
            # and overwrite previously downloaded data
            logging.info('download_url - sleeping %s seconds' % sleep_seconds)
            time.sleep(sleep_seconds)
            download_url(url, local_path, max_retries = max_retries-1, token = token, min_size = min_size)
            return  # success
        else:
            os.remove(local_path)
            raise DownloadError('download_url - failed to download file %s' % url)

    # dump the correct file size to an info file next to the grib file
    # when re-using the GRIB2 file, we check its file size against this record
    # to avoid using partial files
    info_path = local_path + '.size'
    open(ensure_dir(info_path), 'w').write(str(file_size))


def download_many(url_local_pairs, workers=16, aws_only=True, **download_kwargs):
    """
    url_local_pairs: iterable of (url, local_path)
    workers: concurrent downloads
    aws_only: if True, only parallelize s3:// URLs (others can be skipped or handled separately)
    download_kwargs: forwarded to download_url()
    """
    pairs = list(url_local_pairs)

    if aws_only:
        aws_pairs = [(u, p) for (u, p) in pairs if u.startswith("s3://")]
        non_aws_pairs = [(u, p) for (u, p) in pairs if not u.startswith("s3://")]
        if non_aws_pairs:
            logging.info("download_many - %d non-s3 URLs ignored (aws_only=True)", len(non_aws_pairs))
        pairs = aws_pairs

    if not pairs:
        return

    # Threads are appropriate because the heavy work is network + subprocess I/O.
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(download_url, u, p, **download_kwargs): (u, p) for (u, p) in pairs}

        for fut in as_completed(futures):
            u, p = futures[fut]
            try:
                fut.result()
            except Exception as e:
                logging.exception("dowlonad_many - download failed %s -> %s (%s)", u, p, e)