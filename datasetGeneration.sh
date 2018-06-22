
########################################################################
# Generate simulations
########################################################################

mkdir sim

for NUM in {1..20}
do
  rsem-simulated-reads ./rsem-index ./NA12716_7/out.stat/out.model ./NA12716_7/out.isoforms.results 0.0 30000000 ./sim/sim_$NUM --seed $NUM
done
