source ~/.login
/m/linux/matlabR2012b/bin/matlab -nosplash -nodesktop -logfile "~/Collaborations/DIII-D/disruption_database/disruption_database_batch_job.log" -r "d3d_startup; cd('/u/granetzr/Collaborations/DIII-D/disruption_database'); disruption_database(161650:-1:100000); quit"
