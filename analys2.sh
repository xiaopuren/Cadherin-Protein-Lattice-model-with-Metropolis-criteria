for i in {1..10..1}
       do
       cd $i

awk  '{print $1,  $2}' bond.dat  >  homo.dat
awk  '{print $1,  $3}' bond.dat  >  cis.dat
awk  '{print $1,  $4}' bond.dat  >  trans.dat
awk  '{print $1,  $5}' bond.dat  >  cluster.dat
awk  '{print $1,  $6}' bond.dat  >  clusterSize.dat
awk  '{print $1,  $7}' bond.dat  >  maxSize.dat


cd ..
done

paste  1/trans.dat 2/trans.dat 3/trans.dat 4/trans.dat 5/trans.dat 6/trans.dat 7/trans.dat 8/trans.dat 9/trans.dat 10/trans.dat | awk ' { print $1, ($2+ $4+ $6+ $8+ $10+ $12+ $14+ $16+ $18+ $20)/(NF/2)  } ' > ave_trans.dat

paste  1/homo.dat 2/homo.dat 3/homo.dat 4/homo.dat 5/homo.dat 6/homo.dat 7/homo.dat 8/homo.dat 9/homo.dat 10/homo.dat | awk ' { print $1, ($2+ $4+ $6+ $8+ $10+ $12+ $14+ $16+ $18+ $20)/(NF/2)  } ' > ave_homo.dat

paste  1/cis.dat 2/cis.dat 3/cis.dat 4/cis.dat 5/cis.dat 6/cis.dat 7/cis.dat 8/cis.dat 9/cis.dat 10/cis.dat | awk ' { print $1, ($2+ $4+ $6+ $8+ $10+ $12+ $14+ $16+ $18+ $20)/(NF/2)  } ' > ave_cis.dat

paste  1/clusterSize.dat 2/clusterSize.dat 3/clusterSize.dat 4/clusterSize.dat 5/clusterSize.dat 6/clusterSize.dat 7/clusterSize.dat 8/clusterSize.dat 9/clusterSize.dat 10/clusterSize.dat | awk ' { print $1, ($2+ $4+ $6+ $8+ $10+ $12+ $14+ $16+ $18+ $20)/(NF/2)  } ' > ave_clusterSize.dat

paste  1/maxSize.dat 2/maxSize.dat 3/maxSize.dat 4/maxSize.dat 5/maxSize.dat 6/maxSize.dat 7/maxSize.dat 8/maxSize.dat 9/maxSize.dat 10/maxSize.dat | awk ' { print $1, ($2+ $4+ $6+ $8+ $10+ $12+ $14+ $16+ $18+ $20)/(NF/2)  } ' > ave_maxSize.dat
