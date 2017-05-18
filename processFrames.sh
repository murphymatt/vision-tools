##!/bin/bash
mkdir frames
ffmpeg -i videoHD.MOV -r 24 frames/frame_%02d.jpeg
mkdir frames_out
./src/imgpro frames/frame_01.jpeg frames_out/frame_01.jpeg -projectImage frames3/frame_01.jpeg 
ffmpeg -r 24 -i frames_out/frame_%02d.jpeg -c:v libx264 -vf "fps=24,format=yuv420p" out.mp4
