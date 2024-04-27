Run the files with the command 'make [ils/sa]` to generate solutions for the campus data set with output
    Run 'make [ils/sa] Q=1` to have no logging
    Run 'make [ils/sa] L=1` to generate a solution for a dataset with 29 points
    Run 'make [ils/sa] I=X` to change the number of iterations to X (defaults to 50)
    Run 'make [ils/sa] M=X` to generate X amount of solutions and take their averages- will output to sa.csv or ils.csv
        Will actually create a nested loop in the following manner
            1. It will loop from 1 to whatever the iteration limit is (set with the I argument)
            2. The algorithm will run X amount of times with the current iteration limit
            3. The best cost, average cost, and average time for all iteration limits will be written to ils_iteration_limits.csv or sa_iteration_limits.csv

The files `sa.cxx` and `ils.cxx` contain a ton of preprocessor directives/macros that I used to modify the code for 
debugging and data collection without affecting normal execution, but they make the files borderline unreadable in some cases
ils_lean.cxx and sa_lean.cxx contain the source code of the original files, but with most of the directives stripped away. 
I haven't tested these files though, so if they don't run then use the originals