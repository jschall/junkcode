import subprocess
import os
import shutil
from collections import deque
import multiprocessing as mp
import itertools
import threading

imagename = os.path.expanduser('~/livestream/img.jpeg')
streamname = '/dev/video0'

def write_file(imagedata):
    with open(imagename+'.tmp', 'wb') as f:
        f.truncate()
        f.write(imagedata)
    shutil.move(imagename+'.tmp', imagename)

avconv_cmd = ['avconv', '-max_delay', '1000', '-i', streamname, '-max_delay', '1000', '-f', 'image2pipe',  '-r', '1', '-vcodec', 'mjpeg', '-']
avconv = subprocess.Popen(avconv_cmd, stdout = subprocess.PIPE)

#jpeg_header = (0x00, 0x00, 0x00, 0x0C, 0x6A, 0x50, 0x20, 0x20, 0x0D, 0x0A, 0x87, 0x0A) # JPEG-2000
jpeg_header = (0xFF, 0xD8) #JPEG 1998
jpeg_footer = (0xFF, 0xD9)
buf = deque()

worker=None
while True:
    buf.append(ord(avconv.stdout.read(1)))
    while len(buf)>len(jpeg_header):
        buf.popleft()

    if len(buf)>=len(jpeg_header) and tuple(itertools.islice(buf,len(buf)-len(jpeg_header),None)) == jpeg_header:
        # start of image
        while tuple(itertools.islice(buf,len(buf)-len(jpeg_footer),None)) != jpeg_footer:
            buf.append(ord(avconv.stdout.read(1)))
        if worker is None or not worker.isAlive():
            worker = threading.Thread(target=write_file, args=(''.join([chr(b) for b in buf]),))
            worker.start()
        buf.clear()
