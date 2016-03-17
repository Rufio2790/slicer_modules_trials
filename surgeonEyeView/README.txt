Surgeon Eye View Module: this module will be used to verify the correctness of the trajectories implemented. 
Putting the visualization frame on the entry point of the trajectory, I run along the trajectory checking 
slice for slice of the volume. 
    i.	Given two points, define the line between them
    ii.	Create a reference system with the principal axis along the line, and define a plane perpendicular to the line. 
    iii.	Allow the movement slice by slice along the volume, from the entry point to the target point
