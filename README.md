# sTGC.GainDrop

# Run initially like:
python gaindrop.py -f GL3.C3.010_LowCurrent_2019-09-26_16-12.dat.root -g 3
where, 
- the -f argument is the input file
- the -g argument is the gap number
You look at the output 2d plot in order to define the rectangle apparatus and once the rectangle is defined you run again with more arguments

# Then, run like:
python gaindrop.py -f GL3.C3.010_LowCurrent_2019-09-26_16-12.dat.root -g 3 -x 400 -y 700 -c 250,200 -z 200 -l 1
where,
- the -x argument is the width of the rectangle
- the -y argument is the height of the rectangle
- the -c argument is the centre of the rectangle in the form of x0,y0
- the -z argument is the threshold current above which you want to calculate the baseline
The calculation of the area of the blob as well as the average gain-drop in that blob will be done in this step and the rectangle you’ve chosen will be shown on the plot.
