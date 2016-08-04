import subprocess
import os
import shutil
from collections import deque

imagename = os.path.expanduser('~/livestream/img.jpeg')

avconv_cmd = ['avconv', '-i', '/dev/video0', '-f', 'image2pipe',  '-r', '1', '-vcodec', 'mjpeg', '-']
avconv = subprocess.Popen(avconv_cmd, stdout = subprocess.PIPE, stderr = open(os.devnull,'wb'))

#jpeg_header = (0x00, 0x00, 0x00, 0x0C, 0x6A, 0x50, 0x20, 0x20, 0x0D, 0x0A, 0x87, 0x0A) # JPEG-2000
jpeg_header = (0xFF, 0xD8) #JPEG 1998
jpeg_footer = (0xFF, 0xD9)
buf = deque(maxlen=max(len(jpeg_header), len(jpeg_footer)))

while True:
    buf.append(ord(avconv.stdout.read(1)))
    if tuple(buf) == jpeg_header[-len(jpeg_header):]:
        # start of image
        tmp = open(imagename+'.tmp', 'wb')
        tmp.truncate()
        tmp.write(''.join([chr(b) for b in buf]))
        while tuple(buf)[-len(jpeg_footer):] != jpeg_footer:
            buf.append(ord(avconv.stdout.read(1)))
            tmp.write(chr(buf[-1]))
        tmp.flush()
        tmp.close()

        shutil.move(imagename+'.tmp', imagename)
