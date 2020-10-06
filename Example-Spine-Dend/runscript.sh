#!/bin/sh
ugshell -ex reconstructed_spine_wER.lua -grid Spine5_wER.ugx -numRefs 0 -caInflux 0.0086 -tstep 2.5e-6 -endTime 0.015 -setting ryrserca -ryrDensity 0.009 -minDef 1e-12 -outName output -solver GS -vtk -pstep 1e-4

