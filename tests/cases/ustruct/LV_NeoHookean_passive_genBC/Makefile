# Use commands for managing an svFSI simulation and results
clean:
	rm -f -r 4-procs AllData InitialData GenBC.int 
run_svfsi_SimVascular:
	mpirun -np 4 /Users/aaronbrown/Documents/GitHub/svFSI_SimVascular/build/svFSI-build/bin/svFSI svFSI.inp
run_svfsi_vvedula22:
	mpirun -np 4 /Users/aaronbrown/Documents/GitHub/svFSI_vvedula22/build/svFSI-build/bin/svFSI svFSI.inp
run_svfsiplus_SimVascular:
	mpirun -np 4 /Users/aaronbrown/Documents/GitHub/svFSIplus_SimVascular/build/svFSI-build/bin/svFSI svFSI.xml
run_svfsiplus_aabrown100-git:
	mpirun -np 4 /Users/aaronbrown/Documents/GitHub/svFSIplus_aabrown100-git/build/svFSI-build/bin/svFSI svFSI.xml
run_sherlock:
	sbatch ./svFSI_job.sh
stop_sim:
	(cd *-procs/ && touch STOP_SIM)
	