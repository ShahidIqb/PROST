import Model_pred

#path to databases and psiblast
BLAST_DATABASE="$HOME/Downloads/uniref50db/uniref50"
HHBLITS_DATABASE = "$HOME/Downloads/PROST-SEQ_code/data_ddgun/uniclust30_2018_08/uniclust30_2018_08"
BLAST_NUM_THREADS = 12

Model_pred.main(BLAST_DATABASE, BLAST_NUM_THREADS, HHBLITS_DATABASE)
