#!/usr/bin/env python

from ovito.modifiers import HistogramModifier, SelectTypeModifier
from ovito.io import import_file, export_file
import sys

# Get the filename from command-line arguments or use default
filename = sys.argv[1] if len(sys.argv) >= 2 else 'structure.lammpstrj'
pipeline = import_file(filename)
frameNumber = pipeline.source.num_frames
FrameIndexs = [33, 333, 1332, 2331, 3333]

# Apply modifiers to select oxygen atoms and compute histogram
selectO = SelectTypeModifier(property='Particle Type', types={'O'})
histograme = HistogramModifier(bin_count=70, only_selected=True, property='Position.Z', fix_xrange=True, xrange_end=70.3988, xrange_start=0)
pipeline.modifiers.append(selectO)
pipeline.modifiers.append(histograme)

# Export histogram and poscar for specified frames
for frame_index in FrameIndexs:
    print(f"Writing frame {frame_index} to histogram.txt & poscar from {filename}")
    export_file(pipeline, f"{frame_index}.poscar", "vasp", frame=frame_index)
    export_file(pipeline.compute(frame_index), f"{frame_index}histogram.txt", "txt/table", key="histogram[Position.Z]")
