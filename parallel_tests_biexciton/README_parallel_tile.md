Biexciton Script Parallel OMP Style

Notes:
1.) This is only the Rh script. The "qq=" after the max width gamma statement
    in tile_script_top must be updated manually.

2.) The script "tile_set.sh" takes input from "tile_script_top" "*_mid" "*_bottom"
    and user input when the script is executed to create the Fortran file.

    User input when executin ./tile_set.sh:

	a.) Max number: this is the same as the last 'qq' value. Example currently set: 86.
        b.) Number of OMP Threads.

3.) Reasoning behind OMP "Tiling":

	The slow part of the program is when it must search through various criterea and
        sum together different matrix elements. Based on typical output, the output that
	meet the biexcton critera are in the lower ~1/4 of the qq value selected.

	To maximize performance, the tiling script sets up OMP sections for the outermost
	search loop in the biexcton script. One OMP thread is dedicated to the uppder 3/4
	of possible values <qq. The remaining are divided equally across the lower 1/4 block.

	Preliminary tests show this tiling set up is succesful. Work continues to find the 
	best OMP thread amount.

4.) Workflow

   i.) Have all input files needed for serial case and tile_script_* files in the same directory. 
   ii.) Execute ./tile_set.sh. You will need to provide:

		Max number: matches last "qq=" in tile_script_top. For default tile_script_top,
		value is 86.

		Number of OMP Threads.

   iii.) Update the number of OMP threads in the submission script.

   iv.) Submit job.


5.) Bug notes:

	OMP Threads must be greater than 1.
