#!/bin/bash
cd ~/Library/Developer/XCode/DerivedData/HybridFields-hhpmjvywkrgyqlhbtsuaezimjljf/Build/Products/Debug
python z_new.py
mv ./detailedSimulationInfo.txt ./png
cd png
~/../../usr/local/bin/ffmpeg -framerate 6 -start_number 0 -i electricFieldAndParticleAtTime%d.png -vf "scale=720:trunc(ow/a/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p HybrdidFieldsVideo.mp4
