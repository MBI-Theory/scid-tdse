#!/bin/bash
mencoder "mf://*.tiff" -mf fps=6 -o output.avi -ovc lavc \
    -lavcopts vcodec=mpeg4:vbitrate=2400:v4mv:mbd=2:trell:cmp=3:subcmp=3:autoaspect:vpass=1
mencoder "mf://*.tiff" -mf fps=6 -o output.avi -ovc lavc \
    -lavcopts vcodec=mpeg4:vbitrate=2400:v4mv:mbd=2:trell:cmp=3:subcmp=3:autoaspect:vpass=2

