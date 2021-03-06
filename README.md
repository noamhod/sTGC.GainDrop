# sTGC.GainDrop

# Prepare for analysis:
- `mkdir root`
- `mkdir pdf`
- put the root file of the X-ray scan in root/

# Step 1:
`python gaindrop.py -f Your.sTGC.Gap.Xray.Scan.File.root -g 3`

where,
- the -f argument is the input file
- the -g argument is the gap number

Look at the output 2d plot in the pdf/ directory to define the rectangle apparatus and once the rectangle is defined you run again with more arguments as in the instructions below

# Step 2:
`python gaindrop.py -f Your.sTGC.Gap.Xray.Scan.File.root -g 3 -x 400 -y 700 -c 250,200 -z 200 -l 1`

where,
- the -x argument is the width of the rectangle
- the -y argument is the height of the rectangle
- the -c argument is the centre of the rectangle in the form of x0,y0
- the -z argument is the threshold current above which you want to calculate the baseline. Usually it should be a bit smaller than the current value you can read in the corners of the histogram where there's no detector...
- the -l argument is just decifing if the z-axis is plotted in log scale for the 2D plots

The calculation of the area of the blob as well as the average gain-drop in that blob will be done in this step and the rectangle you’ve chosen will be shown on the plots in the pdf/ directory.


# Notes:
- You will probably need to play with the rectangle definition (width, height and center) as well as the threshold currrent to optimise the analysis. This is done by eye in one or two (quick!) iterations.
- This optimisation steps should be repeated untill you see that (a) the blob you see in the plots is reasonably covering the anomaly area, and (b) the entire area of the gap is "blackened" properly, but only the gap area of course (i.e. no significant black area covering the range where there's no detector)
