#!/bin/bash

num_nodi=20				# Numero dei nodi utilizzati
total_sim_time=3600			# Tempo di simulazione
reconfiguration_interval=1000000	# Intervallo dopo cui ho riconfigurazione della rete
lower_bound_speed=1			# Limite minimo della velocità dei nodi
upper_bound_speed=1.86			# Limite massimo della velocità dei nodi
connection_threshold=0.1		# Soglia sopra cui si considera il contatto
x_coord=500				# Lato x dello scenario
y_coord=500				# Lato y dello scenario
rows=20					# Numero di righe della griglia
columns=20				# Numero di colonne della griglia
tx_range=80				# Range trasmissivo dei nodi
rewiring_prob=0.1			# Probabilità di rewiring
groups_number=3				# Numero di gruppi
remaining_prob=0.8			# Probabilita' di rimanere nel gruppo esterno una volta uscito
#seedRNG=9				# Seme del generatore random
travellers_number=1			# Numero dei viaggiatori
travellers_speed=5			# Velocità dei viaggiatori
#travellersProb=0.993
colocationTraces="off"

girmanNewman="off"
deterministic="off"

workingFolder="hcmm"


for seedRNG in "7" #"11" "13" "15" "17" "19" "21" "23"
do
	echo Seed $seedRNG
	
	echo "./hcmm_2_0 -f "$workingFolder$seedRNG" -a $girmanNewman -A $colocationTraces -n $num_nodi -d $deterministic -t $total_sim_time -r $reconfiguration_interval -s $lower_bound_speed -S $upper_bound_speed -p $connection_threshold -X $x_coord -Y $y_coord -R $rows -C $columns -T $tx_range -w $rewiring_prob -G $groups_number -g $seedRNG -c $travellers_number -v $travellers_speed" 

	./hcmm2 -f "$workingFolder$seedRNG" -a $girmanNewman -A $colocationTraces -n $num_nodi -d $deterministic -t $total_sim_time -r $reconfiguration_interval -s $lower_bound_speed -S $upper_bound_speed -p $connection_threshold -X $x_coord -Y $y_coord -R $rows -C $columns -T $tx_range -w $rewiring_prob -G $groups_number -g $seedRNG -c $travellers_number -v $travellers_speed -o $remaining_prob 
	
done
echo ---------------------------------------------------------------

