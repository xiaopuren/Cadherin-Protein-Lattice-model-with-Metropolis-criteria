for i in {1..10..1}
do


    mkdir -p $i

       cd $i

        cp ../main.cpp .
        cp ../run.sh .

sed -i -e 's/#$ -N ad/#$ -N tran3on-'$i'/' run.sh
g++ main.cpp -o out.o
        qsub run.sh

cd ..
done

