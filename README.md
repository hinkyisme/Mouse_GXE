# Mouse_GXE
Gene By Environment Interaction can give clues to personalized diet manage- ment. So to begin, because we measured these mice for multiple phenotypes through 
the experimental timeline, these data are inherently multi-variate. This means that there is likely signal in the data that is nuanced that has interactions 
we can’t see just by doing a simple regression or categorical comparison. Prior to performing any sort of modeling, it would be nice to see what variables 
contribute to the variance seen in the data so we can make informed decisions about what variables (if any) underly the GXE interaction. 
One way to look at what variables(Phenotypes) contribute to this variance is by dimensionality reduction. 
One of the oldest and most used dimensionality reduction methods is a Principle Component Analysis (PCA, click link). 
This method looks to rotate the data to look for orthogonals in the data that show different contributions by different variables (Phenotypes). 
In our case, we will use PCA to reduce the number of dimensions down to just a few, so we can see what the differences are among the variables that contribute to the dataset. 
Therefore, for this experiment we will look if there is a GXE interaction via PCA bi-plots (looking for dimensional variance in Strain but not Environment (Diet)), 
Identify variables that contribute to that variance (e.g., Insulin, Glucose),
Look to see if there are statistical differences in those data by strain to see differential GXE outcomes.
