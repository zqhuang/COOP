for i in `seq 100`
do
    ./LLPW -mask NASYM -isim ${i}
    ./LLPW -mask SASYM -isim ${i}
done    
    
