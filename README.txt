In this folder, there are tools to work with electrical circuits containing diodes. 									
All the scripts have the same coding: the only difference is the electrical model used.									
In DiodeFitting.py, a real diode model is used; in DoubleDiodeFitting.py an additional diode (in opposite direction)	
is added.																												
In this note, there is a brief description of the functionality of these scripts. 										
																														
REQUIREMENTS AND PARAMETERS:
1) This script requires the parameters of the electrical circuit to calculate the IV Curve.
2) The voltage range and the number of acquisition points can be set.
3) The script extracts the values experimentally measured. They are used to:
		a. define the first set of electrical circuit elements (rough, but a starting point);
		b. plot the experimental IV curve
4) It is possible to modify the parameter in order to improve the fitting. A basic knowledge  
   of the impact of the different parameters on the IV curve is therefore necessary.
